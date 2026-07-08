#!/usr/bin/env python3
"""
Batch inference on held-out test set.

For each test tree:
  1. Load root sequence from dataset
  2. Run Algorithm 4 (generate_tree) from that root
  3. Report validity + basic stats vs. ground-truth T1

Usage:
    python scripts/eval_test_set.py \
        --checkpoint checkpoints/best.pt \
        --data       data/train \
        --split      checkpoints/split_indices.json \
        --n-steps    30 \
        --max-trees  20
"""

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import torch
import torch.nn.functional as F
from transformers import AutoTokenizer, EsmForMaskedLM

from src.dataset import TreeDataset
from src.tree_state import TreeState
from src.treeencoder.node_encoder import NodeEncoder
from src.treeencoder.plm_embeddings import ESM2Embedder
from src.treeencoder.structural_features import compute_structural_features
from src.treeencoder.laplacian import compute_laplacian_pe
from src.treeencoder.edges import build_edges
from src.networks import TreeEncoder, RateHeads

AA_VOCAB   = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX  = {aa: i for i, aa in enumerate(AA_VOCAB)}


def load_models(checkpoint, device, max_seq_len):
    ckpt = torch.load(checkpoint, map_location=device, weights_only=False)
    epoch = ckpt.get("epoch", "?")
    print(f"Loaded checkpoint from epoch {epoch}  (val_loss={ckpt.get('val_loss', '?'):.4f})")

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
    L = max_seq_len
    log_rates = torch.zeros(len(sequences), L, 20, dtype=torch.float32, device=device)
    with torch.no_grad():
        tokens = tokenizer(sequences, return_tensors="pt", padding=True,
                           truncation=False).to(device)
        logits = esm_model(**tokens).logits
    seq_lens = tokens["attention_mask"].sum(dim=1)
    for i in range(len(sequences)):
        actual_L = int(seq_lens[i].item()) - 2
        aa_logits = logits[i, 1:actual_L + 1, :][:, aa_token_ids]
        log_probs = F.log_softmax(aa_logits, dim=-1)
        clip = min(actual_L, L)
        log_rates[i, :clip, :] = log_probs[:clip]
    return log_rates


def sequence_identity(a: str, b: str) -> float:
    L = min(len(a), len(b))
    if L == 0:
        return 0.0
    return sum(x == y for x, y in zip(a[:L], b[:L])) / L


def all_valid_aa(seq: str) -> bool:
    return all(c in AA_VOCAB for c in seq)


def generate_one(root_seq, n_steps, max_seq_len, pll_threshold, beta,
                 node_enc, tree_enc, rate_heads, embedder,
                 tokenizer, esm_model, aa_token_ids, device):
    tree = TreeState.root_only(root_seq)
    node_birth_step = {tree.root_id: 0}
    dt = 1.0 / n_steps

    for step in range(n_steps):
        t = step / n_steps
        if not tree.active_leaves:
            break

        node_ids_t  = tree.node_ids
        node_to_idx = {nid: i for i, nid in enumerate(node_ids_t)}
        active_leaves = list(tree.active_leaves)
        active_idx    = [node_to_idx[v] for v in active_leaves]
        node_times_dict = {nid: node_birth_step.get(nid, 0) / n_steps for nid in node_ids_t}

        struct_t = compute_structural_features(tree, node_to_idx).to(device)
        lap_t    = compute_laplacian_pe(tree, node_to_idx, 8, device=device)
        edge_index_t, _, edge_attr_t = build_edges(tree, node_to_idx)
        edge_index_t  = edge_index_t.to(device)
        branch_lens_t = edge_attr_t.squeeze(-1).to(device)

        plm_t = embedder.embed_sequences([tree.node_seqs[nid] for nid in node_ids_t]).to(device)
        with torch.no_grad():
            h_t  = node_enc(plm_t, struct_t, lap_t)
            H_t, _ = tree_enc(h_t, node_ids_t, node_times_dict,
                               edge_index_t, branch_lens_t, t_scalar=t)
            out  = rate_heads(H_t, active_idx)

        active_seqs   = [tree.node_seqs[v] for v in active_leaves]
        log_R0_mut    = get_lm_logits(tokenizer, esm_model, aa_token_ids,
                                      active_seqs, max_seq_len, device)
        log_R_theta   = log_R0_mut + out["mutation_logits"]

        new_node_seqs = dict(tree.node_seqs)

        for i, leaf_id in enumerate(active_leaves):
            seq     = tree.node_seqs[leaf_id]
            seq_len = min(len(seq), max_seq_len)
            new_seq = list(seq)

            for pos in range(seq_len):
                curr_idx = AA_TO_IDX.get(seq[pos], -1)
                if curr_idx < 0:
                    continue
                probs = log_R_theta[i, pos].softmax(-1)
                proposed_idx = torch.multinomial(probs, 1).item()
                if proposed_idx == curr_idx:
                    continue
                delta_pll = (log_R0_mut[i, pos, proposed_idx]
                             - log_R0_mut[i, pos, curr_idx]).item()
                if torch.rand(1).item() < min(1.0, math.exp(beta * delta_pll)) * dt:
                    new_seq[pos] = AA_VOCAB[proposed_idx]
            new_node_seqs[leaf_id] = "".join(new_seq)

            lam  = out["branching_rate"][i].item()
            n_ch = min(int(torch.poisson(torch.tensor(lam * dt)).item()), 2)
            if n_ch > 0:
                child_seqs = [new_node_seqs[leaf_id]] * n_ch
                tree = TreeState(
                    node_ids=tree.node_ids, root_id=tree.root_id,
                    edges=tree.edges, branch_lengths=tree.branch_lengths,
                    node_seqs=new_node_seqs, active_leaves=list(tree.active_leaves),
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

                for child_id in new_children:
                    node_birth_step.setdefault(child_id, step + 1)
                    child_seq = new_node_seqs[child_id]
                    child_len = min(len(child_seq), max_seq_len)
                    aa_idx_c  = torch.tensor(
                        [AA_TO_IDX.get(aa, 20) for aa in child_seq[:child_len]],
                        dtype=torch.long, device=device)
                    valid = aa_idx_c < 20
                    if valid.any():
                        pll = (log_R0_mut[i, :child_len]
                               .gather(-1, aa_idx_c.clamp(max=19).unsqueeze(-1))
                               .squeeze(-1)[valid].mean().item())
                        if pll < pll_threshold:
                            tree = tree.terminate_leaf(child_id)

        tree = TreeState(
            node_ids=tree.node_ids, root_id=tree.root_id,
            edges=tree.edges, branch_lengths=tree.branch_lengths,
            node_seqs=new_node_seqs, active_leaves=list(tree.active_leaves),
        )

    return tree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint",    default="checkpoints/best.pt")
    parser.add_argument("--data",          default="data/train")
    parser.add_argument("--split",         default="checkpoints/split_indices.json")
    parser.add_argument("--n-steps",       type=int,   default=30)
    parser.add_argument("--max-seq-len",   type=int,   default=566)
    parser.add_argument("--pll-threshold", type=float, default=-4.0)
    parser.add_argument("--beta",          type=float, default=1.0)
    parser.add_argument("--max-trees",     type=int,   default=None,
                        help="Cap number of test trees to evaluate (default: all)")
    args = parser.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    node_enc, tree_enc, rate_heads = load_models(args.checkpoint, device, args.max_seq_len)
    embedder = ESM2Embedder(device=device)

    model_id = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_id)
    esm_model = EsmForMaskedLM.from_pretrained(model_id).to(device)
    esm_model.eval()
    for p in esm_model.parameters():
        p.requires_grad = False
    aa_token_ids = torch.tensor(
        [tokenizer.convert_tokens_to_ids(aa) for aa in AA_VOCAB], dtype=torch.long
    )

    dataset = TreeDataset(args.data, max_seq_len=args.max_seq_len)

    with open(args.split) as f:
        split = json.load(f)
    test_indices = split["test"]
    if args.max_trees:
        test_indices = test_indices[:args.max_trees]

    print(f"\nEvaluating {len(test_indices)} test trees  (n_steps={args.n_steps})\n")

    results = []
    for rank, idx in enumerate(test_indices):
        batch = dataset[idx]
        node_ids_T1 = batch["node_ids"]
        seqs_T1     = batch["seqs"]
        root_id_T1  = node_ids_T1[batch["root_index"]]
        root_seq    = seqs_T1[root_id_T1]
        group_num   = batch["group"]

        # ground-truth leaves
        has_children = {p for p, c in batch["edges"]}
        gt_leaves    = [nid for nid in node_ids_T1 if nid not in has_children]
        gt_leaf_seqs = [seqs_T1[nid] for nid in gt_leaves]

        print(f"[{rank+1}/{len(test_indices)}] group={group_num:03d}  "
              f"root_len={len(root_seq)}  gt_nodes={len(node_ids_T1)}  gt_leaves={len(gt_leaves)}")

        try:
            gen_tree = generate_one(
                root_seq, args.n_steps, args.max_seq_len,
                args.pll_threshold, args.beta,
                node_enc, tree_enc, rate_heads, embedder,
                tokenizer, esm_model, aa_token_ids, device,
            )
        except Exception as e:
            print(f"  ERROR: {e}")
            results.append({"group": group_num, "error": str(e)})
            continue

        gen_node_ids = gen_tree.node_ids
        gen_has_ch   = {p for p, c in gen_tree.edges}
        gen_leaves   = [nid for nid in gen_node_ids if nid not in gen_has_ch]
        gen_leaf_seqs = [gen_tree.node_seqs[nid] for nid in gen_leaves]

        # Validity checks
        non_trivial   = len(gen_leaves) >= 2
        valid_seqs    = all(all_valid_aa(s) for s in gen_leaf_seqs)
        bls           = list(gen_tree.branch_lengths.values())
        bls_positive  = all(bl > 0 for bl in bls) if bls else False

        # Best-match sequence identity: each gt leaf matched to closest gen leaf
        if gt_leaf_seqs and gen_leaf_seqs:
            idents = []
            for gs in gt_leaf_seqs:
                best = max(sequence_identity(gs, ps) for ps in gen_leaf_seqs)
                idents.append(best)
            mean_id = sum(idents) / len(idents)
        else:
            mean_id = 0.0

        res = {
            "group":        group_num,
            "gt_nodes":     len(node_ids_T1),
            "gt_leaves":    len(gt_leaves),
            "gen_nodes":    len(gen_node_ids),
            "gen_leaves":   len(gen_leaves),
            "non_trivial":  non_trivial,
            "valid_seqs":   valid_seqs,
            "bls_positive": bls_positive,
            "mean_seq_id":  mean_id,
        }
        results.append(res)

        print(f"  gen_nodes={len(gen_node_ids)}  gen_leaves={len(gen_leaves)}  "
              f"non_trivial={non_trivial}  valid_seqs={valid_seqs}  "
              f"mean_seq_id={mean_id:.3f}")

    # Summary
    ok = [r for r in results if "error" not in r]
    print(f"\n{'='*60}")
    print(f"Results: {len(ok)}/{len(results)} trees generated without error")
    if ok:
        print(f"  Non-trivial (>=2 leaves): {sum(r['non_trivial'] for r in ok)}/{len(ok)}")
        print(f"  Valid AA sequences:       {sum(r['valid_seqs'] for r in ok)}/{len(ok)}")
        print(f"  Positive branch lengths:  {sum(r['bls_positive'] for r in ok)}/{len(ok)}")
        mean_ids = [r['mean_seq_id'] for r in ok]
        print(f"  Mean seq identity:        {sum(mean_ids)/len(mean_ids):.3f}")
        gen_nodes = [r['gen_nodes'] for r in ok]
        print(f"  Gen nodes: min={min(gen_nodes)}  mean={sum(gen_nodes)/len(gen_nodes):.1f}  max={max(gen_nodes)}")

    # Save results
    out_path = Path("checkpoints/eval_test_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
