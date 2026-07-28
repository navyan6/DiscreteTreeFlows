#!/usr/bin/env python3
"""
Aggregate held-out eval: generate a tree from each test root, then report
(1) mutation recovery / conserved retention averaged over roots, and
(2) EVEscape enrichment -- whether the model's mutations land on high-escape
    RBD sites more than a random RBD mutation would, and how that compares to
    the REAL descendants' mutations.

Reuses eval_single_tree's tested generate_tree / load_models / positional_recovery
(so generation is identical to the single-tree smoke test, including the saved
col_entropy and --mutation-rate-scale). EVEscape is position-agnostic, so it
survives generation stochasticity in a way single-tree positional recovery does
not. EVEscape matrix is COVID-spike-RBD only; omit --evescape for H1N1 (recovery
still runs).

Usage:
    python scripts/eval_evescape_enrichment.py \
        --checkpoint checkpoints/covid_v2_entropy/best.pt --data data/covid/test \
        --max-seq-len 1280 --evescape data/covid/evescape_spike_rbd.pt \
        --mutation-rate-scale 0.3 --n-steps 100 --max-trees 20
"""

import argparse
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import torch
from transformers import AutoTokenizer, EsmForMaskedLM

from src.dataset import TreeDataset
from src.tree_state import TreeState
from src.treeencoder.plm_embeddings import ESM2Embedder
from scripts.eval_single_tree import (
    load_models, generate_tree, positional_recovery, get_leaves, seq_identity,
    AA_VOCAB,
)

AA_TO_IDX = {a: i for i, a in enumerate(AA_VOCAB)}


def _mean(xs):
    xs = [x for x in xs if x == x]  # drop nan
    return sum(xs) / len(xs) if xs else float("nan")


def mutation_evescape(root_seq: str, leaf_seq: str, evescape: torch.Tensor, L: int):
    """Per mutation (root->leaf, leaf_aa != root_aa), look up EVEscape[pos, leaf_aa].
    Only positions with a nonzero score (the RBD-scored region) count.
    Returns (list_of_scores, n_muts_total)."""
    scores = []
    n_muts = 0
    for p in range(min(len(root_seq), len(leaf_seq), L)):
        r, g = root_seq[p], leaf_seq[p]
        if g != r and g in AA_TO_IDX:
            n_muts += 1
            s = evescape[p, AA_TO_IDX[g]].item()
            if s != 0.0:
                scores.append(s)
    return scores, n_muts


def gt_leaves_of(batch: dict) -> list[str]:
    parents = {p for p, _ in batch["edges"]}
    return [nid for nid in batch["node_ids"] if nid not in parents]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--max-seq-len", type=int, default=1280)
    ap.add_argument("--evescape", default=None,
                    help="[L,20] EVEscape .pt matrix; enables enrichment (COVID spike only)")
    ap.add_argument("--n-steps", type=int, default=100)
    ap.add_argument("--max-leaves", type=int, default=300)
    ap.add_argument("--branch-rate-scale", type=float, default=6.0)
    ap.add_argument("--mutation-rate-scale", type=float, default=0.3)
    ap.add_argument("--max-trees", type=int, default=20)
    ap.add_argument("--gt-leaves-sampled", type=int, default=30,
                    help="GT leaves per tree to match against gen for recovery")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")
    node_enc, tree_enc, rate_heads, col_entropy = load_models(
        args.checkpoint, device, args.max_seq_len)
    embedder = ESM2Embedder(device=device)
    mid = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(mid)
    esm_model = EsmForMaskedLM.from_pretrained(mid).to(device).eval()
    for p in esm_model.parameters():
        p.requires_grad = False
    aa_token_ids = torch.tensor(
        [tokenizer.convert_tokens_to_ids(a) for a in AA_VOCAB], dtype=torch.long)

    evescape = base_ev = None
    if args.evescape:
        evescape = torch.load(args.evescape, map_location="cpu")
        nz = evescape[evescape != 0]
        base_ev = nz.mean().item()  # mean EVEscape of a random RBD substitution
        print(f"EVEscape {tuple(evescape.shape)}  nonzero={nz.numel()}  "
              f"random-baseline(mean nonzero)={base_ev:.4f}")

    ds = TreeDataset(args.data, max_seq_len=args.max_seq_len)
    n = min(args.max_trees, len(ds)) if args.max_trees else len(ds)
    print(f"Evaluating {n} test trees  (mut_rate={args.mutation_rate_scale} n_steps={args.n_steps})\n")

    recs, rets, idents = [], [], []
    model_ev, gt_ev = [], []            # per-mutation EVEscape scores, pooled
    model_frac_rbd, gt_frac_rbd = [], []  # fraction of mutations that fall in RBD
    per_tree = []
    rng = random.Random(args.seed)

    for i in range(n):
        batch = ds[i]
        root_id = batch["node_ids"][batch["root_index"]]
        root_seq = batch["seqs"][root_id]
        try:
            random.seed(args.seed + i)
            torch.manual_seed(args.seed + i)
            gen = generate_tree(
                root_seq, args.n_steps, args.max_seq_len, args.branch_rate_scale,
                args.max_leaves, args.mutation_rate_scale, node_enc, tree_enc, rate_heads,
                embedder, tokenizer, esm_model, aa_token_ids, device, col_entropy=col_entropy)
        except Exception as e:
            print(f"[{i+1}/{n}] ERROR: {e}")
            continue

        gen_leaves = get_leaves(gen)
        gen_seqs = [gen.node_seqs[g] for g in gen_leaves]
        gt_leaves = gt_leaves_of(batch)
        gt_sample = rng.sample(gt_leaves, min(args.gt_leaves_sampled, len(gt_leaves)))

        # recovery/retention + identity: best-match gen leaf per sampled GT leaf
        t_rec, t_ret, t_id = [], [], []
        for gl in gt_sample:
            gt_seq = batch["seqs"][gl]
            best, best_id = None, -1.0
            for gs in gen_seqs:
                idv = seq_identity(gt_seq, gs)
                if idv > best_id:
                    best_id, best = idv, gs
            if best is None:
                continue
            r = positional_recovery(root_seq, gt_seq, best)
            t_rec.append(r["mut_recovery"]); t_ret.append(r["cons_retention"]); t_id.append(best_id)
        recs.append(_mean(t_rec)); rets.append(_mean(t_ret)); idents.append(_mean(t_id))

        # EVEscape: pool per-mutation scores across gen leaves and (real) GT leaves
        tree_model_ev, tree_gt_ev = [], []
        if evescape is not None:
            for gs in gen_seqs:
                sc, nm = mutation_evescape(root_seq, gs, evescape, args.max_seq_len)
                tree_model_ev += sc
                if nm:
                    model_frac_rbd.append(len(sc) / nm)
            for gl in gt_sample:
                sc, nm = mutation_evescape(root_seq, batch["seqs"][gl], evescape, args.max_seq_len)
                tree_gt_ev += sc
                if nm:
                    gt_frac_rbd.append(len(sc) / nm)
            model_ev += tree_model_ev
            gt_ev += tree_gt_ev

        per_tree.append({"tree": i, "gen_leaves": len(gen_leaves),
                         "mut_recovery": _mean(t_rec), "cons_retention": _mean(t_ret),
                         "identity": _mean(t_id),
                         "model_evescape": _mean(tree_model_ev) if evescape is not None else None,
                         "gt_evescape": _mean(tree_gt_ev) if evescape is not None else None})
        print(f"[{i+1}/{n}] gen_leaves={len(gen_leaves):3d}  "
              f"recovery={_mean(t_rec):.4f}  retention={_mean(t_ret):.4f}  identity={_mean(t_id):.4f}"
              + (f"  model_EVEscape={_mean(tree_model_ev):.4f}" if evescape is not None else ""))

    # ── summary
    print("\n" + "=" * 64)
    print(f"AGGREGATE ({len(recs)} trees)  checkpoint={args.checkpoint}")
    print("=" * 64)
    print(f"  mutation recovery : {_mean(recs):.4f}")
    print(f"  conserved retention: {_mean(rets):.4f}")
    print(f"  best-match identity: {_mean(idents):.4f}")
    summary = {"checkpoint": args.checkpoint, "n_trees": len(recs),
               "mut_recovery": _mean(recs), "cons_retention": _mean(rets),
               "identity": _mean(idents)}
    if evescape is not None:
        m_ev, g_ev = _mean(model_ev), _mean(gt_ev)
        print(f"\n  EVEscape enrichment (higher = mutations at more escape-prone RBD sites):")
        print(f"    model mutations   : mean EVEscape = {m_ev:.4f}   ({len(model_ev)} scored muts, "
              f"{_mean(model_frac_rbd)*100:.1f}% of muts in RBD)")
        print(f"    real GT mutations : mean EVEscape = {g_ev:.4f}   ({len(gt_ev)} scored muts, "
              f"{_mean(gt_frac_rbd)*100:.1f}% of muts in RBD)")
        print(f"    random RBD baseline: mean EVEscape = {base_ev:.4f}")
        print(f"    model - random    : {m_ev - base_ev:+.4f}   (model vs GT: {m_ev - g_ev:+.4f})")
        summary.update({"model_evescape": m_ev, "gt_evescape": g_ev,
                        "random_baseline_evescape": base_ev,
                        "model_frac_rbd": _mean(model_frac_rbd),
                        "gt_frac_rbd": _mean(gt_frac_rbd)})

    out = args.out or f"checkpoints/eval_enrichment_{Path(args.checkpoint).parent.name}.json"
    Path(out).write_text(json.dumps({"summary": summary, "per_tree": per_tree}, indent=2))
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main()
