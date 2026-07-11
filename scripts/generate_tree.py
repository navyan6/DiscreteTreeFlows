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
    N, L = len(sequences), max_seq_len
    log_rates = torch.zeros(N, L, 20, dtype=torch.float32, device=device)
    with torch.no_grad():
        tokens = tokenizer(sequences, return_tensors="pt", padding=True,
                           truncation=False).to(device)
        logits = esm_model(**tokens).logits          # [N, max_len+2, vocab]
    seq_lens = tokens["attention_mask"].sum(dim=1)
    for i in range(N):
        actual_L = int(seq_lens[i].item()) - 2
        aa_logits = logits[i, 1:actual_L + 1, :][:, aa_token_ids]
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
    seq = args.root_seq.upper()
    nt_chars = set("ACGTU")
    nt_frac = sum(c in nt_chars for c in seq) / max(len(seq), 1)
    if nt_frac > 0.85:
        raise ValueError(
            f"root-seq looks like a nucleotide sequence ({nt_frac:.0%} ACGTU). "
            "Translate to amino acids first."
        )
    aa_frac = sum(c in AA_VOCAB for c in seq) / max(len(seq), 1)
    print(f"Root sequence: {len(seq)} aa  ({aa_frac:.0%} standard AA)")

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

        # R0 log-rates computed first so mutation head can condition on them
        active_seqs = [tree.node_seqs[v] for v in active_leaves]
        log_R0_mut  = get_lm_logits(tokenizer, esm_model, aa_token_ids,
                                     active_seqs, args.max_seq_len, device)

        # NodeEncoder → TreeEncoder → RateHeads
        with torch.no_grad():
            h_t  = node_enc(plm_t, struct_t, lap_t)
            H_t, _ = tree_enc(h_t, node_ids_t, node_times_dict,
                               edge_index_t, branch_lens_t, t_scalar=t)
            out  = rate_heads(H_t, active_idx, log_R0_mut)
        # out["log_R_theta_mut"] = log_R0 + c_θ  [n_active, L, 20]

        # Sample events for each active leaf
        new_node_seqs = dict(tree.node_seqs)

        for i, leaf_id in enumerate(active_leaves):
            seq     = tree.node_seqs[leaf_id]
            seq_len = min(len(seq), args.max_seq_len)

            # ── Mutations (M-H: model-guided proposal + ESM fitness gate) ──
            new_seq = list(seq)
            for pos in range(seq_len):
                curr_idx = AA_TO_IDX.get(seq[pos], -1)
                if curr_idx < 0:
                    continue
                probs = out["log_R_theta_mut"][i, pos].softmax(-1)
                probs_mut = probs.clone()
                probs_mut[curr_idx] = 0.0
                total = probs_mut.sum().item()
                if total > 0 and torch.rand(1).item() < total * dt:
                    new_seq[pos] = AA_VOCAB[torch.multinomial(probs_mut / total, 1).item()]
            new_node_seqs[leaf_id] = "".join(new_seq)

            # ── Branch ──
            lam  = out["branching_rate"][i].item() * args.branch_rate_scale
            p_branch = 1.0 - math.exp(-max(0.0, lam) * dt)
            n_ch = 2 if torch.rand(1).item() < p_branch else 0
            if n_ch > 0:
                child_seqs = [new_node_seqs[leaf_id]] * n_ch
                tree = TreeState(
                    node_ids=tree.node_ids, root_id=tree.root_id,
                    edges=tree.edges, branch_lengths=tree.branch_lengths,
                    node_seqs=new_node_seqs,
                    active_leaves=list(tree.active_leaves),
                )
                tree = tree.branch_node(leaf_id, child_seqs)

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

                # ── ESM fitness gate: terminate new children with very low PLL ──
                for child_id in new_children:
                    node_birth_step.setdefault(child_id, step + 1)
                    child_seq = new_node_seqs[child_id]
                    child_len = min(len(child_seq), args.max_seq_len)
                    aa_idx_c = torch.tensor(
                        [AA_TO_IDX.get(aa, 20) for aa in child_seq[:child_len]],
                        dtype=torch.long, device=device,
                    )
                    valid = aa_idx_c < 20
                    if valid.any():
                        child_pll = (
                            log_R0_mut[i, :child_len]
                            .gather(-1, aa_idx_c.clamp(max=19).unsqueeze(-1))
                            .squeeze(-1)[valid]
                            .mean()
                            .item()
                        )
                        if child_pll < args.pll_threshold:
                            tree = tree.terminate_leaf(child_id)

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

    # Save all node sequences as FASTA
    fasta_path = out_path.with_suffix(".fasta")
    has_children = {p for p, c in tree.edges}
    with open(fasta_path, "w") as f:
        for nid in tree.node_ids:
            tag = "root" if nid == tree.root_id else ("leaf" if nid not in has_children else "internal")
            f.write(f">{nid}|{tag}\n{tree.node_seqs[nid]}\n")
    print(f"Sequences saved to {fasta_path}")

    return tree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint",    required=True,
                        help="Path to trained checkpoint (e.g. checkpoints/best.pt)")
    parser.add_argument("--root-seq",      required=True,
                        help="Root amino acid sequence")
    parser.add_argument("--n-steps",       type=int,   default=50)
    parser.add_argument("--output",        default="generated_tree.nwk")
    parser.add_argument("--max-seq-len",   type=int,   default=566)
    parser.add_argument("--pll-threshold",    type=float, default=-100.0,
                        help="Terminate new child if ESM PLL < this (nats/position); -100 disables gate")
    parser.add_argument("--beta",             type=float, default=1.0,
                        help="MH acceptance temperature (higher = stricter ESM fitness gate)")
    parser.add_argument("--branch-rate-scale", type=float, default=6.0,
                        help="Multiply model branching rate by this at inference (corrects lam≈1 → lam≈6)")
    args = parser.parse_args()
    generate_tree(args)


if __name__ == "__main__":
    main()
