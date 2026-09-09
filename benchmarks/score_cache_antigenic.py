#!/usr/bin/env python3
"""Score lit/antigenic hotspot mut fraction from coverage tree caches (C.3 / D.2).

Antigenic = mean over gen leaves of hotspot_mut_frac(root, leaf), same definition
as eval_evescape_enrichment.lit_hotspot_mut_frac. Optional aa|hit if --test-data.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.eval_evescape_enrichment import (  # noqa: E402
    hotspot_mut_frac,
    load_hotspot_mask,
)
from scripts.eval_single_tree import get_leaves, positional_recovery  # noqa: E402
from src.tree_state import TreeState  # noqa: E402


VIRUS_MASK = {
    "covid": ("results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt", 1280),
    "h1n1": ("results/flu_mutfreq_vs_lit/mut_hotspot_mask_h1_nmicrobiol_lit.pt", 566),
    "h3n2": ("results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt", 566),
    "hiv_geo": ("results/hiv_env_mask/mut_hotspot_mask_hiv_env_v1v5.pt", 900),
    "hiv_temporal": ("results/hiv_env_mask/mut_hotspot_mask_hiv_env_v1v5.pt", 900),
}


def load_tree(path: Path) -> TreeState:
    return TreeState.from_dict(json.loads(path.read_text()))


def score_gen_antigenic(root_dir: Path, mask: torch.Tensor, L: int, k_max: int) -> dict:
    trees = []
    for i in range(k_max):
        p = root_dir / f"tree_{i:04d}.json"
        if not p.exists():
            break
        trees.append(load_tree(p))
    if not trees:
        return {}
    root_id = trees[0].root_id
    root_seq = trees[0].node_seqs.get(root_id, "")
    ant = []
    n_leaves = 0
    for tree in trees:
        for lab in get_leaves(tree):
            seq = tree.node_seqs.get(lab)
            if not seq:
                continue
            n_leaves += 1
            if seq == root_seq:
                continue
            f = hotspot_mut_frac(root_seq, seq, mask, L)
            if f == f:
                ant.append(f)
    return {
        "root_id": root_dir.name,
        "n_trees": len(trees),
        "n_gen_leaves": n_leaves,
        "lit_hotspot_mut_frac": mean(ant) if ant else float("nan"),
        "n_mut_leaves_scored": len(ant),
    }


def score_aa_hit(
    root_dir: Path,
    gt_leaf_seqs: list[str],
    k_max: int,
) -> dict:
    from benchmarks.metrics import sequences as S

    trees = []
    for i in range(k_max):
        p = root_dir / f"tree_{i:04d}.json"
        if not p.exists():
            break
        trees.append(load_tree(p))
    if not trees or not gt_leaf_seqs:
        return {}
    root_id = trees[0].root_id
    root_seq = trees[0].node_seqs.get(root_id, "")
    gen_pool = []
    for tree in trees:
        for lab in get_leaves(tree):
            seq = tree.node_seqs.get(lab)
            if seq:
                gen_pool.append(seq)
    aa_hits, mut_rec, cons_ret = [], [], []
    for g in gt_leaf_seqs:
        best = max(gen_pool, key=lambda x: S.identity(g, x))
        pr = positional_recovery(root_seq, g, best)
        if pr["aa_acc_given_hit"] == pr["aa_acc_given_hit"]:
            aa_hits.append(pr["aa_acc_given_hit"])
        if pr["mut_recovery"] == pr["mut_recovery"]:
            mut_rec.append(pr["mut_recovery"])
        if pr["cons_retention"] == pr["cons_retention"]:
            cons_ret.append(pr["cons_retention"])
    return {
        "aa_acc_given_hit": mean(aa_hits) if aa_hits else float("nan"),
        "mut_recovery": mean(mut_rec) if mut_rec else float("nan"),
        "cons_retention": mean(cons_ret) if cons_ret else float("nan"),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--virus", required=True, choices=sorted(VIRUS_MASK))
    ap.add_argument("--cache-dir", required=True,
                    help="coverage_cache_* directory (contains METHOD/N16/ROOT)")
    ap.add_argument("--methods", nargs="+", default=None)
    ap.add_argument("--N", type=int, default=16)
    ap.add_argument("--K-max", type=int, default=100)
    ap.add_argument("--per-root-csv", default=None,
                    help="C.3 *_per_root.csv to map root_id→horizon_bucket")
    ap.add_argument("--test-data", default=None,
                    help="If set, also score aa|hit / mut / cons vs GT leaves")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    mask_path, L = VIRUS_MASK[args.virus]
    mask, meta = load_hotspot_mask(str(ROOT / mask_path), L)
    if mask is None:
        raise SystemExit(f"missing mask {mask_path}")

    cache = Path(args.cache_dir)
    if not cache.is_absolute():
        cache = ROOT / cache

    horizon_by_root: dict[str, str] = {}
    if args.per_root_csv:
        prp = Path(args.per_root_csv)
        if not prp.is_absolute():
            prp = ROOT / prp
        with open(prp) as f:
            for row in csv.DictReader(f):
                horizon_by_root[row["root_id"]] = row["horizon_bucket"]

    gt_by_root: dict[str, list[str]] = {}
    if args.test_data:
        from benchmarks.heldout.build_examples import build_examples
        from benchmarks.metrics import trees as T
        from benchmarks.run_table import rebuild_target

        for ex in build_examples(ROOT / args.test_data, args.N, seed=0):
            tgt = rebuild_target(ex)
            leaves = [tgt.node_seqs[lab] for lab in T.leaf_labels(tgt)
                      if tgt.node_seqs.get(lab)]
            gt_by_root[ex["root_id"]] = leaves

    methods = args.methods
    if methods is None:
        methods = sorted(d.name for d in cache.iterdir() if d.is_dir())

    per_root = []
    for method in methods:
        method_dir = cache / method / f"N{args.N}"
        if not method_dir.is_dir():
            print(f"WARN: missing {method_dir}")
            continue
        for rd in sorted(d for d in method_dir.iterdir() if d.is_dir()):
            rec = score_gen_antigenic(rd, mask, L, args.K_max)
            if not rec:
                continue
            rec["method"] = method
            rec["horizon_bucket"] = horizon_by_root.get(rd.name)
            if args.test_data and rd.name in gt_by_root:
                rec.update(score_aa_hit(rd, gt_by_root[rd.name], args.K_max))
            per_root.append(rec)
            print(json.dumps({k: rec[k] for k in rec if k != "n_gen_leaves"}))

    def agg(rows: list[dict], key: str) -> float | None:
        vals = [r[key] for r in rows if key in r and r[key] == r[key]]
        return mean(vals) if vals else None

    by_mh: dict[tuple, list] = defaultdict(list)
    for r in per_root:
        by_mh[(r.get("horizon_bucket") or "all", r["method"])].append(r)

    cells = []
    for (h, m), rows in sorted(by_mh.items()):
        cells.append({
            "horizon_bucket": h,
            "method": m,
            "n_roots": len(rows),
            "lit_hotspot_mut_frac": agg(rows, "lit_hotspot_mut_frac"),
            "aa_acc_given_hit": agg(rows, "aa_acc_given_hit"),
            "mut_recovery": agg(rows, "mut_recovery"),
            "cons_retention": agg(rows, "cons_retention"),
        })

    summary = {
        "virus": args.virus,
        "cache": str(cache),
        "mask": mask_path,
        "K_max": args.K_max,
        "n_roots_scored": len(per_root),
        "cells": cells,
        "per_root": per_root,
    }
    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))
    print("WROTE", out)
    print(json.dumps({"cells": cells}, indent=2))


if __name__ == "__main__":
    main()
