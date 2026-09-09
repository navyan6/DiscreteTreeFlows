#!/usr/bin/env python3
"""
Appendix C.2 Ab metrics from Rodriguez 82 samples.

Fills / refreshes for pLM prior + treesbm_ab_oas:
  - Cov@K at absolute Hamming e=2 (default K∈{100,500})
  - Viral-style clade_recall + mean_min_edit at K=100 (same defs as coverage_curves)

Mut. recall in C.2 Ab rows stays CDR mut. recall from Track C (not recomputed here).

Usage:
  python scripts/eval_ab_c2_from_samples.py \\
    --models plm_prior,treesbm_ab_oas \\
    --K-list 100,500 \\
    --out benchmarks/results/tables/table_c2_ab_metrics.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from benchmarks.coverage_curves import (  # noqa: E402
    clade_recall_at_k,
    emerging_clade_fingerprints,
)
from benchmarks.metrics import sequences as S  # noqa: E402
from scripts.ab_t6_metrics import coverage_at_k_abs  # noqa: E402
from scripts.eval_ab_t6_from_samples import (  # noqa: E402
    _gen_leaf_aa,
    _leaf_aa,
    _load_generated,
    _load_trees,
    _mean,
)
from src.tree_state import TreeState  # noqa: E402


def _rod_to_tree_state(tree: dict) -> TreeState | None:
    aa = tree.get("true_aa_sequences") or {}
    if not aa:
        return None
    bl: dict[tuple[str, str], float] = {}
    for k, v in (tree.get("branch_lengths") or {}).items():
        if "->" in str(k):
            a, b = str(k).split("->", 1)
            bl[(a, b)] = float(v)
        else:
            return None
    edges = [tuple(e) for e in (tree.get("edges") or [])]
    node_ids = tree.get("node_ids") or sorted({n for e in edges for n in e})
    try:
        return TreeState(
            node_ids=list(node_ids),
            root_id=tree["root_id"],
            edges=edges,  # type: ignore[arg-type]
            branch_lengths=bl,
            node_seqs={nid: aa[nid] for nid in node_ids if nid in aa},
        )
    except Exception:
        return None


def _pool_gen_leaves(tree: dict, rollouts: list[dict]) -> list[str]:
    out: list[str] = []
    for item in rollouts:
        seqs = item.get("sequences") or {}
        out.extend(
            _gen_leaf_aa(
                tree,
                seqs,
                gen_root_id=item.get("root_id"),
                free_topology=bool(item.get("free_topology", False)),
            )
        )
    return out


def evaluate_model(
    trees: list[dict],
    gens: dict[str, list[dict]],
    model: str,
    *,
    k_list: list[int],
    viral_k: int,
    seed: int,
) -> dict:
    cov_lists = {k: [] for k in k_list}
    clade_list: list[float] = []
    edit_list: list[float] = []
    pool_sizes: list[int] = []
    n_eligible: list[int] = []
    per_family = []

    for tree in trees:
        fid = tree["family_id"]
        root, gt_leaves = _leaf_aa(tree)
        rollouts = gens.get(fid, [])
        if not root or not gt_leaves or not rollouts:
            continue
        gen_leaves = _pool_gen_leaves(tree, rollouts)
        if not gen_leaves:
            continue
        pool_sizes.append(len(gen_leaves))
        fam: dict = {
            "family_id": fid,
            "n_gt_leaves": len(gt_leaves),
            "n_gen_leaves": len(gen_leaves),
            "n_rollouts": len(rollouts),
        }
        for k in k_list:
            c = coverage_at_k_abs(
                gt_leaves, gen_leaves, e=2, k=k, seed=seed + hash(fid) % 10_000
            )
            fam[f"coverage_at_{k}_e2"] = c
            cov_lists[k].append(c)

        # Viral-style clade / min dist at viral_k (match C.2 K=100 convention)
        pool = list(gen_leaves)
        if len(pool) > viral_k:
            import random

            rng = random.Random(seed + hash(fid) % 10_000)
            rng.shuffle(pool)
            pool = pool[:viral_k]
        edits = [S.min_hamming(g, pool) for g in gt_leaves]
        mean_edit = sum(edits) / len(edits) if edits else float("nan")
        fam["mean_min_edit"] = mean_edit
        edit_list.append(mean_edit)

        gt_state = _rod_to_tree_state(tree)
        if gt_state is not None:
            fps = emerging_clade_fingerprints(gt_state, root)
            cr, n_el = clade_recall_at_k(fps, pool, root)
            fam["clade_recall"] = cr
            fam["n_eligible_clades"] = n_el
            n_eligible.append(n_el)
            if cr == cr:
                clade_list.append(cr)
        else:
            fam["clade_recall"] = float("nan")
            fam["n_eligible_clades"] = 0
        per_family.append(fam)

    summary = {
        "model": model,
        "n_families": len(per_family),
        "viral_K": viral_k,
        "median_n_gen_leaves": (
            sorted(pool_sizes)[len(pool_sizes) // 2] if pool_sizes else float("nan")
        ),
        "n_families_pool_lt_500": sum(1 for n in pool_sizes if n < 500),
        "clade_recall": _mean(clade_list),
        "mean_min_edit": _mean(edit_list),
        "mean_n_eligible_clades": _mean([float(x) for x in n_eligible]),
    }
    for k in k_list:
        summary[f"coverage_at_{k}_e2"] = _mean(cov_lists[k])
    return {"summary": summary, "per_family": per_family}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", type=Path, default=ROOT / "antibody_benchmark/data/processed/benchmark_trees.jsonl")
    ap.add_argument("--samples-dir", type=Path, default=ROOT / "antibody_benchmark/results/samples")
    ap.add_argument("--models", default="plm_prior,treesbm_ab_oas")
    ap.add_argument("--K-list", default="100,500")
    ap.add_argument("--viral-K", type=int, default=100, help="K for clade_recall / mean_min_edit")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--md-out", type=Path, default=None)
    args = ap.parse_args()

    k_list = [int(x) for x in args.K_list.split(",") if x.strip()]
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    trees = _load_trees(args.trees)
    print(f"trees={len(trees)} models={models} K={k_list} viral_K={args.viral_K}", flush=True)

    results = {}
    rows = []
    for model in models:
        gens = _load_generated(args.samples_dir, model)
        print(f"model={model} families_with_samples={len(gens)}", flush=True)
        if not gens:
            results[model] = {"summary": {"model": model, "n_families": 0, "status": "no_samples"}, "per_family": []}
            continue
        block = evaluate_model(
            trees, gens, model, k_list=k_list, viral_k=args.viral_K, seed=args.seed
        )
        results[model] = block
        rows.append(block["summary"])
        s = block["summary"]
        print(
            f"  n={s['n_families']} cov100={s.get('coverage_at_100_e2')} "
            f"cov500={s.get('coverage_at_500_e2')} clade={s.get('clade_recall')} "
            f"edit={s.get('mean_min_edit')} pool_med={s.get('median_n_gen_leaves')} "
            f"pool_lt500={s.get('n_families_pool_lt_500')}",
            flush=True,
        )

    payload = {
        "protocol": (
            "C.2 Ab: Cov@K e2 from pooled gen leaf AAs; clade_recall + mean_min_edit "
            "via coverage_curves defs on Rod.82 GT TreeState (AA). Mut = CDR (Track C)."
        ),
        "trees": str(args.trees),
        "samples_dir": str(args.samples_dir),
        "K_list": k_list,
        "viral_K": args.viral_K,
        "models": results,
        "table_rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n")

    md = args.md_out or args.out.with_suffix(".md")
    cols = [
        "model",
        "n_families",
        *[f"coverage_at_{k}_e2" for k in k_list],
        "clade_recall",
        "mean_min_edit",
        "median_n_gen_leaves",
        "n_families_pool_lt_500",
    ]

    def fmt(v):
        if v is None:
            return "—"
        if isinstance(v, float):
            if isinstance(v, float) and math.isnan(v):
                return "NaN"
            return f"{v:.4f}"
        return str(v)

    lines = [
        "# Appendix C.2 — Ab Cov@500 / clade / min dist (Rod.82)",
        "",
        f"Source: `{args.out}`",
        "",
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for r in rows:
        lines.append("| " + " | ".join(fmt(r.get(c)) for c in cols) + " |")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- **Cov@K e2**: absolute Hamming ≤2; pool = up to K gen leaf AAs across rollouts.",
            "- **Clade recall / Min dist**: viral `coverage_curves` defs at K=100 on Rod.82 GT trees.",
            "- **Mut. recall** for C.2 Ab rows remains **CDR mut. recall** from Track C.",
            "- If `n_families_pool_lt_500` > 0, Cov@500 uses the available pool (<500) for those families.",
            "",
        ]
    )
    md.write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.out}")
    print(f"Wrote {md}")


if __name__ == "__main__":
    main()
