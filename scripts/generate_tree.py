#!/usr/bin/env python3
"""
Algorithm 4: Controlled tree generation using trained TreeSBM model.

Starting from a root sequence, generates a bifurcating phylogenetic tree by
sampling from R_theta(T, T', t) = R0(T, T') * exp(c_theta(T, t)).

Mutation sampling:
  For each position pos, the probability of a mutation in step dt is
  (1 - p_current_aa) * dt, where p_current_aa = softmax(log_R_theta_mut)[pos, current_aa].
  Positions the model is confident about mutate rarely; uncertain ones mutate more.

Bifurcating constraint: sampled child count is clamped to max 2.

Usage:
    python scripts/generate_tree.py \\
        --checkpoint checkpoints/best.pt \\
        --root-seq ACDEFG... \\
        --n-steps 50 \\
        --output generated_tree.nwk
"""

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import torch
import torch.nn.functional as F
from transformers import AutoTokenizer, EsmForMaskedLM

from src.tree_state import TreeState
from src.treeencoder.node_encoder import NodeEncoder
from src.treeencoder.plm_embeddings import ESM2Embedder
from src.treeencoder.structural_features import compute_structural_features
from src.treeencoder.laplacian import compute_laplacian_pe
from src.treeencoder.edges import build_edges
from src.networks import TreeEncoder, RateHeads

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_VOCAB)}


def load_checkpoint(path, device, max_seq_len=566):
    ckpt = torch.load(path, map_location=device, weights_only=False)
    node_enc  = NodeEncoder(d_plm=320, d_struct=3, d_laplacian=8, d_node=128).to(device)
    tree_enc  = TreeEncoder(d_model=128, n_layers=4, n_heads=8, dropout=0.1).to(device)
    r_heads   = RateHeads(d_model=128, max_seq_len=max_seq_len).to(device)
    node_enc.load_state_dict(ckpt["node_enc"])
    tree_enc.load_state_dict(ckpt["tree_enc"])
    r_heads.load_state_dict(ckpt["rate_heads"])
    node_enc.eval(); tree_enc.eval(); r_heads.eval()
    for m in [node_enc, tree_enc, r_heads]:
        for p in m.parameters():
            p.requires_grad = False
    return node_enc, tree_enc, r_heads


def get_lm_logits(tokenizer, esm_model, aa_token_ids, sequences, max_seq_len, device):
#get ESM logits
    N, L = len(sequences), max_seq_len
    log_rates = torch.zeros(N, L, 20, dtype=torch.float32, device=device)
    for i, seq in enumerate(sequences):
        with torch.no_grad():
            tokens = tokenizer(seq, return_tensors="pt").to(device)
            logits = esm_model(**tokens).logits          # [1, actual_L+2, 33]
        actual_L = int(tokens["attention_mask"].sum().item()) - 2
        aa_logits = logits[0, 1:actual_L + 1, :][:, aa_token_ids]  # [actual_L, 20]
        log_probs = F.log_softmax(aa_logits, dim=-1)
        clip = min(actual_L, L)
        log_rates[i, :clip, :] = log_probs[:clip]
    return log_rates


def tree_to_newick(tree: TreeState) -> str:
#treestate to newich str
    children_map: dict[str, list[str]] = {}
    for p, c in tree.edges:
        children_map.setdefault(p, []).append(c)

    def _fmt(nid: str) -> str:
        bl = tree.branch_lengths.get(
            next(((p, nid) for p, c in tree.edges if c == nid), (None, None)),
            0.0
        )
        children = children_map.get(nid, [])
        if not children:
            return f"{nid}:{bl:.6f}"
        inner = ",".join(_fmt(c) for c in children)
        return f"({inner}){nid}:{bl:.6f}"

    return _fmt(tree.root_id) + ";"


def generate_tree(args):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    node_enc, tree_enc, rate_heads = load_checkpoint(args.checkpoint, device, args.max_seq_len)

    # ESM-2-8M 
    embedder = ESM2Embedder(device=device)

    # ESM-2-8M 
    model_id = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_id)
    esm_model = EsmForMaskedLM.from_pretrained(model_id).to(device)
    esm_model.eval()
    for p in esm_model.parameters():
        p.requires_grad = False
    aa_token_ids = torch.tensor(
        [tokenizer.convert_tokens_to_ids(aa) for aa in AA_VOCAB], dtype=torch.long
    )

    # t=0: root-only tree
    tree = TreeState.root_only(args.root_seq)
    # Track birth step for each node (used as proxy for calendar time in causal mask)
    node_birth_step: dict[str, int] = {tree.root_id: 0}

    n_steps = args.n_steps
    dt = 1.0 / n_steps

    for step in range(n_steps):
        t = step / n_steps

        if not tree.active_leaves:
            print(f"Step {step}: no active leaves, stopping early")
            break

        node_ids_t    = tree.node_ids
        node_to_idx   = {nid: i for i, nid in enumerate(node_ids_t)}
        active_leaves = list(tree.active_leaves)
        active_idx    = [node_to_idx[v] for v in active_leaves]

        node_times_dict = {
            nid: node_birth_step.get(nid, 0) / n_steps for nid in node_ids_t
        }

        # Tree features for T_t
        struct_t = compute_structural_features(tree, node_to_idx).to(device)
        lap_t    = compute_laplacian_pe(tree, node_to_idx, 8, device=device)
        edge_index_t, _, edge_attr_t = build_edges(tree, node_to_idx)
        edge_index_t  = edge_index_t.to(device)
        branch_lens_t = edge_attr_t.squeeze(-1).to(device)

        # PLM embeddings [N, 320] for current sequences
        seqs_list = [tree.node_seqs[nid] for nid in node_ids_t]
        plm_t = embedder.embed_sequences(seqs_list).to(device)

        # NodeEncoder → TreeEncoder → RateHeads
        with torch.no_grad():
            h_t  = node_enc(plm_t, struct_t, lap_t)
            H_t, _ = tree_enc(h_t, node_ids_t, node_times_dict,
                               edge_index_t, branch_lens_t, t_scalar=t)
            out  = rate_heads(H_t, active_idx)

        # R0 log-rates for active leaf sequences
        active_seqs = [tree.node_seqs[v] for v in active_leaves]
        log_R0_mut   = get_lm_logits(tokenizer, esm_model, aa_token_ids,
                                      active_seqs, args.max_seq_len, device)
        log_R_theta_mut = log_R0_mut + out["mutation_logits"]   # [n_active, L, 20]

        # Sample events for each active leaf
        new_node_seqs = dict(tree.node_seqs)

        for i, leaf_id in enumerate(active_leaves):
            seq     = tree.node_seqs[leaf_id]
            seq_len = min(len(seq), args.max_seq_len)

            #  PLL check 
            aa_idx_i = torch.tensor(
                [AA_TO_IDX.get(aa, 20) for aa in seq[:seq_len]],
                dtype=torch.long, device=device,
            )
            valid = aa_idx_i < 20
            if valid.any():
                pll_i = (
                    log_R0_mut[i, :seq_len]
                    .gather(-1, aa_idx_i.clamp(max=19).unsqueeze(-1))
                    .squeeze(-1)[valid]
                    .mean()
                    .item()
                )
                if pll_i < args.pll_threshold:
                    tree = tree.terminate_leaf(leaf_id)
                    continue

            #  Stop 
            if torch.bernoulli(out["stop_prob"][i]).item():
                tree = tree.terminate_leaf(leaf_id)
                continue

            #  Mutations (Metropolis-Hastings: model-guided proposal + ESM fitness gate) ──
            new_seq = list(seq)
            for pos in range(seq_len):
                curr_idx = AA_TO_IDX.get(seq[pos], -1)
                if curr_idx < 0:
                    continue

                # Propose from model distribution (log R0 + c_theta)
                probs = log_R_theta_mut[i, pos].softmax(-1)
                proposed_idx = torch.multinomial(probs, 1).item()
                if proposed_idx == curr_idx:
                    continue

                # Accept/reject by ESM fitness (Metropolis-Hastings)
                delta_pll = (log_R0_mut[i, pos, proposed_idx]
                             - log_R0_mut[i, pos, curr_idx]).item()
                accept_prob = min(1.0, math.exp(args.beta * delta_pll))

                if torch.rand(1).item() < accept_prob * dt:
                    new_seq[pos] = AA_VOCAB[proposed_idx]
            new_node_seqs[leaf_id] = "".join(new_seq)

            # Branch 
            # Clamp to 2: phylogenetic trees are bifurcating
            lam  = out["branching_rate"][i].item()
            n_ch = min(int(torch.poisson(torch.tensor(lam * dt)).item()), 2)
            if n_ch > 0:
                child_seqs = [new_node_seqs[leaf_id]] * n_ch
                # Rebuild tree with updated sequences before branching
                tree = TreeState(
                    node_ids=tree.node_ids, root_id=tree.root_id,
                    edges=tree.edges, branch_lengths=tree.branch_lengths,
                    node_seqs=new_node_seqs,
                    active_leaves=list(tree.active_leaves),
                )
                tree = tree.branch_node(leaf_id, child_seqs)

                # Assign predicted branch length to new edges
                bl_pred = out["branch_length"][i].item()
                new_children = tree.get_children(leaf_id)
                new_bls = {(leaf_id, c): bl_pred for c in new_children}
                tree = TreeState(
                    node_ids=tree.node_ids, root_id=tree.root_id,
                    edges=tree.edges,
                    branch_lengths={**tree.branch_lengths, **new_bls},
                    node_seqs=tree.node_seqs,
                    active_leaves=list(tree.active_leaves),
                )
                new_node_seqs = dict(tree.node_seqs)

                for child_id in new_children:
                    node_birth_step.setdefault(child_id, step + 1)

        # Flush sequence updates for non-branching leaves
        tree = TreeState(
            node_ids=tree.node_ids, root_id=tree.root_id,
            edges=tree.edges, branch_lengths=tree.branch_lengths,
            node_seqs=new_node_seqs,
            active_leaves=list(tree.active_leaves),
        )

        print(
            f"Step {step + 1:03d}/{n_steps}  "
            f"nodes={len(tree.node_ids)}  "
            f"active_leaves={len(tree.active_leaves)}"
        )

    print(f"\nFinal tree: {len(tree.node_ids)} nodes, {len(tree.active_leaves)} leaves")

    nwk = tree_to_newick(tree)
    out_path = Path(args.output)
    out_path.write_text(nwk)
    print(f"Saved to {out_path}")
    return tree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint",    required=True,
                        help="Path to trained checkpoint (e.g. checkpoints/best.pt)")
    parser.add_argument("--root-seq",      required=True,
                        help="Root amino acid sequence")
    parser.add_argument("--n-steps",       type=int,   default=50)
    parser.add_argument("--output",        default="generated_tree.nwk")
    parser.add_argument("--pll-threshold", type=float, default=-2.5,
                        help="Terminate branch if mean ESM PLL < this (nats/position)")
    parser.add_argument("--max-seq-len",   type=int,   default=566)
    parser.add_argument("--beta",          type=float, default=1.0,
                        help="MH acceptance temperature (higher = stricter ESM fitness gate)")
    args = parser.parse_args()
    generate_tree(args)


if __name__ == "__main__":
    main()
