#!/usr/bin/env python3
"""
Track A1 — full future generation, benchmarked on real held-out trees.

For each test group: build the GT TreeState, sample K future trees from its root,
and score the generated set against the single real tree with the metric library
(sequence recovery + tree-shape agreement + sequence-topology coupling). Since
real data gives one GT tree per root, distributional distances are reported as
gen-vs-GT shape ratios rather than distribution-vs-distribution.

Runs where a checkpoint + ESM are available (cluster/GPU):

    python benchmarks/track_a.py \
        --checkpoint checkpoints/h3n2_temporal/best.pt \
        --data data/h3n2/test --K 10 --n-groups 20 \
        --out benchmarks/results/track_a_h3n2_temporal.json
"""

import argparse
import json
import random
import sys
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from src.dataset import TreeDataset
from benchmarks.generation import TreeSBMGenerator, gt_treestate
from benchmarks.metrics import trees as T
from benchmarks.metrics import sequences as S
from benchmarks.metrics import matching as MT


def _safe_ratio(a: float, b: float) -> float:
    return a / b if b else float("nan")


def score_group(gt, root_seq, gen_trees, max_gt=60, max_gen_pool=800, seed=0):
    rng = random.Random(seed)
    gt_leaves = T.leaf_labels(gt)
    gt_leaf_seqs = [gt.node_seqs[l] for l in gt_leaves]

    # pool generated leaves across the K trees
    gen_pool = [t.node_seqs[l] for t in gen_trees for l in T.leaf_labels(t)]
    gen_pool_sub = (rng.sample(gen_pool, max_gen_pool)
                    if len(gen_pool) > max_gen_pool else gen_pool)
    gt_sub = rng.sample(gt_leaf_seqs, max_gt) if len(gt_leaf_seqs) > max_gt else gt_leaf_seqs

    # ── sequence recovery (blind: leaves not identity-matched)
    prf_sub = S.mutation_pr_f1(gen_pool, gt_leaf_seqs, root_seq, level="substitution")
    prf_site = S.mutation_pr_f1(gen_pool, gt_leaf_seqs, root_seq, level="site")
    best_of_k = mean(S.best_of_k_identity(g, gen_pool_sub) for g in gt_sub)
    cov = S.coverage_at_k(gt_sub, gen_pool_sub, eps_frac=0.02)

    # positional recovery: each GT leaf -> its nearest generated leaf, GT root anchor
    mut_rec, cons_ret = [], []
    for g in gt_sub:
        best = max(gen_pool_sub, key=lambda x: S.identity(g, x))
        pr = S.positional_recovery(root_seq, g, best)
        if pr["mut_total"] or pr["cons_total"]:
            mut_rec.append(pr["mut_recovery"])
            cons_ret.append(pr["cons_retention"])

    # ── tree-shape agreement (gen mean vs GT)
    gt_shape = {"sackin": T.sackin_index(gt), "colless": T.colless_index(gt),
                "cherries": T.cherry_count(gt), "height": T.topological_height(gt),
                "mean_bl": mean(gt.branch_lengths.values()) if gt.branch_lengths else 0.0}
    g_sackin = mean(T.sackin_index(t) for t in gen_trees)
    g_colless = mean(T.colless_index(t) for t in gen_trees)
    g_cherry = mean(T.cherry_count(t) for t in gen_trees)
    g_height = mean(T.topological_height(t) for t in gen_trees)
    g_bl = mean(mean(t.branch_lengths.values()) for t in gen_trees if t.branch_lengths)
    coupling = mean(MT.seq_patristic_correlation(t) for t in gen_trees)

    def _fin(v):
        return [x for x in v if x == x]

    return {
        "mut_recovery": mean(_fin(mut_rec)) if _fin(mut_rec) else float("nan"),
        "cons_retention": mean(_fin(cons_ret)) if _fin(cons_ret) else float("nan"),
        "best_of_k_identity": best_of_k,
        "coverage@k_eps2pct": cov,
        "mut_precision": prf_sub["precision"], "mut_recall": prf_sub["recall"],
        "mut_f1": prf_sub["f1"], "site_recall": prf_site["recall"],
        "sackin_ratio": _safe_ratio(g_sackin, gt_shape["sackin"]),
        "colless_ratio": _safe_ratio(g_colless, gt_shape["colless"]),
        "cherry_ratio": _safe_ratio(g_cherry, gt_shape["cherries"]),
        "height_ratio": _safe_ratio(g_height, gt_shape["height"]),
        "branch_len_ratio": _safe_ratio(g_bl, gt_shape["mean_bl"]),
        "seq_patristic_corr": coupling,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--data", required=True, help="dir of processed test groups")
    ap.add_argument("--groups", nargs="+", type=int, default=None)
    ap.add_argument("--n-groups", type=int, default=20, help="if --groups not given, use the first N")
    ap.add_argument("--K", type=int, default=10)
    ap.add_argument("--n-steps", type=int, default=50)
    ap.add_argument("--max-leaves", type=int, default=400)
    ap.add_argument("--branch-rate-scale", type=float, default=6.0)
    ap.add_argument("--mutation-rate-scale", type=float, default=0.04)
    ap.add_argument("--max-seq-len", type=int, default=566)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    dataset = TreeDataset(args.data, max_seq_len=args.max_seq_len)
    avail = sorted(dataset.groups)
    groups = args.groups or avail[: args.n_groups]
    print(f"{len(avail)} groups available; evaluating {len(groups)}: {groups}")

    gen = TreeSBMGenerator(args.checkpoint, max_seq_len=args.max_seq_len)

    per_group = {}
    for g in groups:
        idx = next((i for i in range(len(dataset)) if dataset.groups[i] == g), None)
        if idx is None:
            print(f"[{g}] not found, skipping"); continue
        gt, root_seq = gt_treestate(dataset[idx])
        gen_trees = gen.generate_k(
            root_seq, args.K, n_steps=args.n_steps, max_leaves=args.max_leaves,
            branch_rate_scale=args.branch_rate_scale,
            mutation_rate_scale=args.mutation_rate_scale, base_seed=args.seed,
        )
        res = score_group(gt, root_seq, gen_trees, seed=args.seed)
        per_group[g] = res
        print(f"[{g}] recovery={res['mut_recovery']:.3f} cons={res['cons_retention']:.3f} "
              f"cov@k={res['coverage@k_eps2pct']:.3f} mutF1={res['mut_f1']:.3f} "
              f"BLratio={res['branch_len_ratio']:.2f} couple={res['seq_patristic_corr']:.2f}")

    # aggregate (mean over groups, ignoring NaN)
    keys = next(iter(per_group.values())).keys() if per_group else []
    agg = {k: mean([v[k] for v in per_group.values() if v[k] == v[k]]) for k in keys}
    print("\n=== Track A1 aggregate (mean over groups) ===")
    for k, v in agg.items():
        print(f"  {k:22s} {v:.4f}")

    out = {"checkpoint": args.checkpoint, "data": args.data, "K": args.K,
           "groups": list(per_group.keys()), "aggregate": agg, "per_group": per_group}
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(out, indent=2))
        print(f"\nSaved {args.out}")


if __name__ == "__main__":
    main()
