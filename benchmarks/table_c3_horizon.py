#!/usr/bin/env python3
"""
Table C.3 — coverage stratified by genetic horizon (train H terciles).

Locked protocol:
  - Horizon buckets = short / med / long from train-set ``h_buckets`` (1/3, 2/3).
  - True conditioning: ``generate(root_seq, N, H)`` with each example's own GT
    mean root-to-tip H (not post-hoc filtering of generated trees).
  - Stratify metrics by which train tercile that example's H falls in.
  - Same roots as Table 5 (COVID / H1N1 / HIV) or Table 4 (H3N2); N=16.
  - Methods: NeutralBD, ARTreeFormer-adapted, TreeSBM (virus ckpt).
  - Metrics: Coverage@K e-abs e1/e2/e3/e5, mut_recovery, clade_recall,
    mean_min_edit. Tree-KL omitted (not cheap at K=100).

Examples:

  python benchmarks/table_c3_horizon.py \\
      --virus covid --pin-groups 4,5,7,9,10 \\
      --checkpoint checkpoints/covid_v5_mutrec/best.pt --max-seq-len 1280 \\
      --methods neutral_bd artreeformer_adapted treesbm \\
      --K-list 10 100 --out benchmarks/results/table_c3_covid.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from benchmarks.coverage_curves import (  # noqa: E402
    DEFAULT_E_LIST,
    _leaf_seqs,
    _parse_e_list,
    build_baseline_methods,
    emerging_clade_fingerprints,
    generate_or_load_trees,
    score_pool,
)
from benchmarks.heldout.build_examples import (  # noqa: E402
    build_examples,
    h_buckets,
    list_groups,
)
from benchmarks.metrics import sequences as S  # noqa: E402
from benchmarks.run_table import rebuild_target  # noqa: E402


# Locked Table 5 / Table 4 root groups (1-indexed group ids).
DEFAULT_PIN_GROUPS = {
    "covid": [4, 5, 7, 9, 10],          # Brazil diverse (job 7575015)
    "h1n1": [1, 4, 5, 6, 7],            # geo roots (job 7575016)
    "h3n2": [2],                        # Table 4 group_002 (5 roots from one tree)
    "hiv_geo": None,                    # diverse one-per-group (same filter as T5)
    "hiv_temporal": None,
}

VIRUS_DEFAULTS = {
    "covid": {
        "test_data": "data/covid/test",
        "train_data": "data/covid/train",
        "params": "benchmarks/results/params_covid.json",
        "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
        "max_seq_len": 1280,
        "max_roots_per_tree": 1,
    },
    "h1n1": {
        "test_data": "data/h1n1/test",
        "train_data": "data/h1n1/train",
        "params": "benchmarks/results/params_h1n1.json",
        "checkpoint": "checkpoints/h1n1_v2_lit_hotspot/best.pt",
        "max_seq_len": 566,
        "max_roots_per_tree": 1,
    },
    "h3n2": {
        "test_data": "data/h3n2/test",
        "train_data": "data/h3n2/train",
        "params": "benchmarks/results/params.json",
        "checkpoint": "checkpoints/h3n2_v3_lit_hotspot/best.pt",
        "max_seq_len": 566,
        "max_roots_per_tree": 5,  # Table 4: 5 roots from group_002
    },
    "hiv_geo": {
        "test_data": "data/hiv_geo/test",
        "train_data": "data/hiv_geo/train",
        "params": "benchmarks/results/params_hiv_geo.json",
        "checkpoint": "checkpoints/hiv_geo_v1/best.pt",
        "max_seq_len": 900,
        "max_roots_per_tree": 1,
    },
    "hiv_temporal": {
        "test_data": "data/hiv_temporal/test",
        "train_data": "data/hiv_temporal/train",
        "params": "benchmarks/results/params_hiv_temporal.json",
        "checkpoint": "checkpoints/hiv_temporal_v1/best.pt",
        "max_seq_len": 900,
        "max_roots_per_tree": 1,
    },
}


def horizon_label(H: float, lo: float, hi: float) -> str:
    if H <= lo:
        return "short"
    if H <= hi:
        return "med"
    return "long"


def collect_train_H(train_dir: Path, N: int, seed: int,
                    max_roots_per_tree: int = 1) -> list[float]:
    """Train-set mean RTT values used for tercile cuts (same builder as eval)."""
    raw = build_examples(
        train_dir, N, seed=seed, max_roots_per_tree=max_roots_per_tree)
    return [float(ex["H"]) for ex in raw]

def select_examples(
    test_dir: Path,
    N: int,
    seed: int,
    max_roots: int,
    max_roots_per_tree: int,
    pin_groups: list[int] | None,
    diverse_filter: bool = True,
) -> list[dict]:
    raw = build_examples(
        test_dir, N, seed=seed, max_roots_per_tree=max_roots_per_tree)

    if pin_groups:
        want = set(int(g) for g in pin_groups)
        pinned = [ex for ex in raw if int(ex["group"]) in want]
        # Preserve pin order, then root order within group
        order = {g: i for i, g in enumerate(pin_groups)}
        pinned.sort(key=lambda e: (order.get(int(e["group"]), 999), e["root_id"]))
        if len(pinned) > max_roots:
            pinned = pinned[:max_roots]
        if not pinned:
            raise SystemExit(
                f"no examples for pin-groups={pin_groups} in {test_dir}")
        return pinned

    examples: list[dict] = []
    skipped = 0
    for ex in raw:
        leaves = list(ex["target_seqs"].values())
        if not leaves:
            continue
        if diverse_filter:
            mean_edit = mean(S.hamming(ex["root_seq"], s) for s in leaves)
            L = max(len(ex["root_seq"]), 1)
            if mean_edit <= 0 or mean_edit > 0.05 * L:
                skipped += 1
                continue
        examples.append(ex)
        if len(examples) >= max_roots:
            break
    if len(examples) < max_roots:
        print(f"WARN: only {len(examples)} diverse roots "
              f"(skipped={skipped}); falling back to first {max_roots}")
        examples = raw[:max_roots]
    return examples


def aggregate_rows(rows: list[dict], keys: list[str]) -> dict:
    out: dict = {}
    for k in keys:
        vals = [r[k] for r in rows if r.get(k) == r.get(k)]
        out[k] = mean(vals) if vals else float("nan")
    return out


def run(args) -> None:
    virus = args.virus
    defaults = VIRUS_DEFAULTS[virus]
    test_dir = ROOT / (args.test_data or defaults["test_data"])
    train_dir = ROOT / (args.train_data or defaults["train_data"])
    params_path = ROOT / (args.params or defaults["params"])
    ckpt = args.checkpoint or defaults["checkpoint"]
    max_seq_len = args.max_seq_len or defaults["max_seq_len"]
    max_roots_per_tree = (
        args.max_roots_per_tree
        if args.max_roots_per_tree is not None
        else defaults["max_roots_per_tree"]
    )

    if not list_groups(test_dir):
        raise SystemExit(f"no processed groups in {test_dir}")
    if not params_path.exists():
        raise SystemExit(
            f"missing {params_path}; fit with: "
            f"python benchmarks/fit_params.py --train-data {train_dir} "
            f"--out {params_path}")

    params = json.loads(params_path.read_text())
    N = args.N
    ks = list(args.K_list)
    e_list = _parse_e_list(args.e_list) if args.e_list else [1, 2, 3, 5]

    # --- train H terciles ---
    train_H = collect_train_H(train_dir, N, seed=args.seed, max_roots_per_tree=1)
    if len(train_H) < 3:
        raise SystemExit(f"too few train H values ({len(train_H)}) for terciles")
    lo, hi = h_buckets(train_H)
    cuts = {
        "virus": virus,
        "N": N,
        "n_train_H": len(train_H),
        "H_min": min(train_H),
        "H_max": max(train_H),
        "H_mean": mean(train_H),
        "tercile_lo": lo,
        "tercile_hi": hi,
        "buckets": {
            "short": f"H <= {lo:.6g}",
            "med": f"{lo:.6g} < H <= {hi:.6g}",
            "long": f"H > {hi:.6g}",
        },
        "protocol": (
            "generate(root, N, H=example_mean_RTT); "
            "stratify by train h_buckets terciles"
        ),
    }
    cuts_path = ROOT / args.cuts_out
    cuts_path.parent.mkdir(parents=True, exist_ok=True)
    cuts_path.write_text(json.dumps(cuts, indent=2))
    print(f"train H terciles: lo={lo:.6g} hi={hi:.6g} "
          f"(n={len(train_H)}) -> {cuts_path}")

    # --- select roots ---
    pin = args.pin_groups
    if pin is None and virus in DEFAULT_PIN_GROUPS:
        pin = DEFAULT_PIN_GROUPS[virus]
    if pin is not None and isinstance(pin, str):
        pin = [int(x) for x in pin.split(",") if x.strip()]
    examples = select_examples(
        test_dir, N, args.seed, args.max_roots, max_roots_per_tree,
        pin_groups=pin, diverse_filter=not args.no_diverse_filter,
    )
    for ex in examples:
        ex["horizon_bucket"] = horizon_label(float(ex["H"]), lo, hi)

    root_meta = [
        {
            "group": e["group"],
            "root_id": e["root_id"],
            "H": e["H"],
            "horizon_bucket": e["horizon_bucket"],
            "mean_edit": mean(
                S.hamming(e["root_seq"], s) for s in e["target_seqs"].values()
            ) if e["target_seqs"] else None,
        }
        for e in examples
    ]
    print(f"roots ({len(examples)}):")
    for r in root_meta:
        print(f"  group={r['group']:03d} root={r['root_id']} "
              f"H={r['H']:.4g} bucket={r['horizon_bucket']} "
              f"mean_edit={r['mean_edit']:.2f}")

    # --- methods ---
    class _A:
        pass
    margs = _A()
    margs.methods = args.methods
    margs.N = [N]
    try:
        margs.train_data = str(train_dir.relative_to(ROOT))
    except ValueError:
        margs.train_data = str(train_dir)
    margs.checkpoint = ckpt
    margs.n_steps = args.n_steps
    margs.max_seq_len = max_seq_len
    margs.seq_model = args.seq_model
    margs.r0_backend = None
    margs.r0_model = None
    margs.fitness_beta = None
    margs.ablate_bridge = False
    margs.ablate_tree_context = False
    margs.ablate_branch_length_head = False
    margs.ablate_internal_node_seqs = False
    margs.ablate_site_entropy = False
    margs.branching_mode = "learned"
    margs.ref_lambda = 1.0
    methods = build_baseline_methods(margs, params, esm=None)
    if not methods:
        raise SystemExit("no methods available")
    print(f"methods={[m.name for m in methods]} Ks={ks} e_list={e_list}")

    # generation args object reused by coverage_curves helpers
    cache_dir = Path(args.cache_dir) if args.cache_dir else None
    if cache_dir and not cache_dir.is_absolute():
        cache_dir = ROOT / cache_dir
    rescore_from = Path(args.rescore_from) if args.rescore_from else None
    if rescore_from and not rescore_from.is_absolute():
        rescore_from = ROOT / rescore_from

    gargs = _A()
    gargs.K_max = max(ks)
    gargs.seed = args.seed
    gargs.cache_dir = str(cache_dir) if cache_dir else None
    gargs.rescore_from = str(rescore_from) if rescore_from else None
    gargs.eps_frac = args.eps_frac
    gargs.max_gt_leaves = args.max_gt_leaves
    gargs.max_gen_pool = args.max_gen_pool
    gargs.min_clade_size = 2
    gargs.min_shared_muts = 1
    gargs.e_list = e_list

    # --- generate + score ---
    metric_keys = [
        "coverage", "mean_min_edit", "mean_min_edit_frac",
        "site_recall", "mut_f1", "mut_recovery", "cons_retention",
        "aa_acc_given_hit", "best_of_k_identity",
        "clade_recall", "n_eligible_clades",
    ]
    for e in e_list:
        metric_keys += [
            f"coverage_obs_e{e}",
            f"coverage_obs_unique_e{e}",
            f"frac_gen_e{e}",
        ]

    per_root_fields = [
        "virus", "method", "N", "K", "horizon_bucket",
        "group", "root_id", "H", "tercile_lo", "tercile_hi",
        *metric_keys, "n_trees", "runtime_gen_sec",
    ]
    agg_fields = [
        "virus", "method", "N", "K", "horizon_bucket",
        "n_roots", "n_trees_scored", "tercile_lo", "tercile_hi",
        *metric_keys, "runtime_gen_sec",
    ]

    per_root_rows: list[dict] = []
    # (method, K, bucket) -> list of metric dicts
    buckets: dict[tuple, list[dict]] = {}
    gen_secs_by_method: dict[str, float] = {}

    for method in methods:
        gen_secs_by_method[method.name] = 0.0
        for ex in examples:
            target = rebuild_target(ex)
            gt_leaves = _leaf_seqs(target)
            if len(gt_leaves) != N:
                print(f"  skip {ex['root_id']}: gt leaves={len(gt_leaves)} != N={N}")
                continue
            bucket = ex["horizon_bucket"]
            trees, secs = generate_or_load_trees(
                method, ex, gargs, N,
                cache_dir=cache_dir,
                rescore_from=rescore_from,
            )
            gen_secs_by_method[method.name] += secs
            if not trees:
                continue
            fps = emerging_clade_fingerprints(target, ex["root_seq"])
            pool: list[str] = []
            next_i = 0
            for K in ks:
                while next_i < min(K, len(trees)):
                    pool.extend(_leaf_seqs(trees[next_i]))
                    next_i += 1
                metrics = score_pool(
                    gt_leaves, ex["root_seq"], pool,
                    eps_frac=args.eps_frac,
                    max_gt=args.max_gt_leaves,
                    max_gen=args.max_gen_pool,
                    seed=args.seed + K,
                    e_list=e_list,
                    gt_tree=target,
                    fingerprints=fps,
                )
                row = {
                    "virus": virus,
                    "method": method.name,
                    "N": N,
                    "K": K,
                    "horizon_bucket": bucket,
                    "group": ex["group"],
                    "root_id": ex["root_id"],
                    "H": ex["H"],
                    "tercile_lo": lo,
                    "tercile_hi": hi,
                    "n_trees": len(trees),
                    "runtime_gen_sec": secs,
                    **metrics,
                }
                per_root_rows.append(row)
                key = (method.name, K, bucket)
                buckets.setdefault(key, []).append(metrics)
                cov_e2 = metrics.get("coverage_obs_e2", float("nan"))
                print(
                    f"  [{method.name}] g={ex['group']:03d} root={ex['root_id']} "
                    f"bucket={bucket} K={K} "
                    f"cov@e2={cov_e2 if cov_e2 == cov_e2 else float('nan'):.3f} "
                    f"edit={metrics['mean_min_edit']:.1f} "
                    f"mut_rec={metrics['mut_recovery']:.3f}"
                )

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=agg_fields, extrasaction="ignore")
        w.writeheader()
        for (method_name, K, bucket), rows in sorted(
            buckets.items(), key=lambda x: (x[0][0], x[0][2], x[0][1])
        ):
            agg = aggregate_rows(rows, metric_keys)
            w.writerow({
                "virus": virus,
                "method": method_name,
                "N": N,
                "K": K,
                "horizon_bucket": bucket,
                "n_roots": len(rows),
                "n_trees_scored": K,
                "tercile_lo": lo,
                "tercile_hi": hi,
                "runtime_gen_sec": gen_secs_by_method.get(method_name, 0.0),
                **agg,
            })
    print(f"wrote aggregate {out}")

    per_root_out = ROOT / args.per_root_out
    per_root_out.parent.mkdir(parents=True, exist_ok=True)
    with open(per_root_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=per_root_fields, extrasaction="ignore")
        w.writeheader()
        for row in per_root_rows:
            w.writerow(row)
    print(f"wrote per-root {per_root_out}")

    # also dump root meta + cuts together
    meta_path = out.with_suffix(".meta.json")
    meta_path.write_text(json.dumps({
        "cuts": cuts,
        "roots": root_meta,
        "methods": [m.name for m in methods],
        "K_list": ks,
        "e_list": e_list,
        "checkpoint": ckpt,
        "max_seq_len": max_seq_len,
        "pin_groups": pin,
    }, indent=2))
    print(f"wrote meta {meta_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--virus", required=True, choices=list(VIRUS_DEFAULTS))
    ap.add_argument("--test-data", default=None)
    ap.add_argument("--train-data", default=None)
    ap.add_argument("--params", default=None)
    ap.add_argument("--checkpoint", default=None)
    ap.add_argument("--max-seq-len", type=int, default=None)
    ap.add_argument("--N", type=int, default=16)
    ap.add_argument("--K-list", type=int, nargs="+", default=[10, 100])
    ap.add_argument("--max-roots", type=int, default=5)
    ap.add_argument("--max-roots-per-tree", type=int, default=None)
    ap.add_argument("--pin-groups", default=None,
                    help="Comma-separated 1-indexed groups (default: locked T5/T4)")
    ap.add_argument("--no-diverse-filter", action="store_true")
    ap.add_argument("--methods", nargs="+",
                    default=["neutral_bd", "artreeformer_adapted", "treesbm"],
                    choices=["neutral_bd", "plm_prior", "artreeformer_adapted", "treesbm"])
    ap.add_argument("--e-list", nargs="+", default=["1", "2", "3", "5"])
    ap.add_argument("--eps-frac", type=float, default=0.02)
    ap.add_argument("--n-steps", type=int, default=50)
    ap.add_argument("--seq-model", default="JTT")
    ap.add_argument("--max-gt-leaves", type=int, default=60)
    ap.add_argument("--max-gen-pool", type=int, default=800)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--rescore-from", default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None,
                    help="Aggregate-by-bucket CSV")
    ap.add_argument("--per-root-out", default=None)
    ap.add_argument("--cuts-out", default=None)
    args = ap.parse_args()

    if args.out is None:
        args.out = f"benchmarks/results/table_c3_{args.virus}_N{args.N}.csv"
    if args.per_root_out is None:
        args.per_root_out = (
            f"benchmarks/results/table_c3_{args.virus}_N{args.N}_per_root.csv")
    if args.cuts_out is None:
        args.cuts_out = (
            f"benchmarks/results/tables/table_c3_{args.virus}_h_cuts.json")
    # Namespace expects e_list as list for _parse_e_list
    run(args)


if __name__ == "__main__":
    main()
