#!/usr/bin/env python3
"""
Table 5 antibody maturation eval harness (scaffolding).

Wires:
  - CDR mut recall / SHM load / terminal diversity from scripts/ab_t6_metrics.py
  - Lineage RF via dendropy (when both GT + gen Newick exist)
  - Coverage@100 via benchmarks.metrics.sequences.coverage_at_k (K=100 pool)

Baselines from benchmarks.methods.ab_baselines — Neutral SHM and AR raise
unless fallback flags set. Does not invent paper numbers.

Usage (after data/ab_clones/test exists):
  python scripts/eval_ab_maturation.py \\
    --data data/ab_clones/test \\
    --methods plm_prior \\
    --out results/ab_t5/eval_smoke.json \\
    --max-groups 5 --K 10
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.ab_t6_metrics import (  # noqa: E402
    cdr_mask_imgt_stub,
    cdr_mut_recall,
    coverage_at_k_abs,
    shm_load_error,
    terminal_diversity_error,
)


def load_group(dir_path: Path, gname: str):
    fa = dir_path / f"{gname}.fasta"
    nwk = dir_path / f"{gname}.nwk"
    root_fa = dir_path / f"{gname}_root.fasta"
    if not fa.is_file() or not nwk.is_file():
        return None
    leaves = []
    header = None
    parts: list[str] = []
    with fa.open() as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    leaves.append("".join(parts))
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
        if header is not None:
            leaves.append("".join(parts))
    root = ""
    if root_fa.is_file():
        seqs = []
        p = []
        for line in root_fa.open():
            if line.startswith(">"):
                if p:
                    seqs.append("".join(p))
                    p = []
            else:
                p.append(line.strip())
        if p:
            seqs.append("".join(p))
        root = seqs[0] if seqs else ""
    if not root and leaves:
        # fallback: first leaf (bad; flagged)
        root = leaves[0]
    return {"leaves": leaves, "root": root, "nwk": nwk.read_text().strip(), "root_is_placeholder": not root_fa.is_file()}


def rf_distance(nwk_a: str, nwk_b: str) -> float:
    try:
        import dendropy
        from dendropy.calculate import treecompare
    except ImportError:
        return float("nan")
    tns = dendropy.TaxonNamespace()
    ta = dendropy.Tree.get(data=nwk_a, schema="newick", taxon_namespace=tns)
    tb = dendropy.Tree.get(data=nwk_b, schema="newick", taxon_namespace=tns)
    ta.encode_bipartitions()
    tb.encode_bipartitions()
    return float(treecompare.symmetric_difference(ta, tb))


def coverage_at_100(gt_leaves, gen_leaves) -> float:
    try:
        from benchmarks.metrics import sequences as S
        return float(S.coverage_at_k(gt_leaves, gen_leaves, eps_frac=0.02))
    except Exception:
        return float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--methods", default="plm_prior", help="Comma list: neutral_shm,plm_prior,ar_tree_edit")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max-groups", type=int, default=20)
    ap.add_argument("--K", type=int, default=10, help="Trees / samples per root for Coverage pool")
    ap.add_argument("--N", type=int, default=32)
    ap.add_argument("--H", type=float, default=1.0)
    ap.add_argument("--birth", type=float, default=1.0)
    ap.add_argument("--death", type=float, default=0.5)
    ap.add_argument("--allow-neutral-bd-fallback", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    groups = sorted({p.name.split(".")[0] for p in args.data.glob("group_*.fasta")})
    groups = groups[: args.max_groups]
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    # Lazy baseline construction
    from benchmarks.methods.ab_baselines import build_ab_baseline

    esm_logits = None
    built = {}
    for m in methods:
        if m.startswith("plm"):
            if esm_logits is None:
                try:
                    from benchmarks.run_table import ESM, _cuda
                    esm_logits = ESM("cuda" if _cuda() else "cpu").lm_logits
                except Exception as e:
                    print(f"WARN: ESM unavailable ({e}); plm_prior will fail", flush=True)
            try:
                built[m] = build_ab_baseline(
                    m, args.birth, args.death, esm_logits=esm_logits,
                    allow_neutral_bd_fallback=args.allow_neutral_bd_fallback,
                )
            except Exception as e:
                built[m] = e
        else:
            try:
                built[m] = build_ab_baseline(
                    m, args.birth, args.death,
                    allow_neutral_bd_fallback=args.allow_neutral_bd_fallback,
                )
            except Exception as e:
                built[m] = e

    rows = []
    for gi, gname in enumerate(groups):
        g = load_group(args.data, gname)
        if not g:
            continue
        mask = cdr_mask_imgt_stub(len(g["root"]))
        for m in methods:
            meth = built.get(m)
            rec = {
                "group": gname,
                "method": m,
                "n_leaves_gt": len(g["leaves"]),
                "status": "ok",
            }
            if isinstance(meth, Exception) or meth is None:
                rec["status"] = f"baseline_unavailable: {meth}"
                rows.append(rec)
                continue
            gen_leaves = []
            gen_nwks = []
            try:
                for k in range(args.K):
                    gt = meth.generate(g["root"], min(args.N, len(g["leaves"])), args.H, args.seed + 1000 * gi + k)
                    leaves = [gt.tree.node_seqs[l] for l in gt.tree.active_leaves]
                    gen_leaves.extend(leaves)
                    # topology RF needs newick — skip if TreeState only
                rec["cdr_mut_recall"] = cdr_mut_recall(g["root"], gen_leaves, g["leaves"], mask)
                rec["shm_load_error"] = shm_load_error(g["root"], gen_leaves, g["leaves"])
                rec["terminal_diversity_error"] = terminal_diversity_error(gen_leaves, g["leaves"])
                rec["coverage_at_100"] = coverage_at_k_abs(
                    g["leaves"],
                    gen_leaves,
                    e=2,
                    k=min(100, max(1, len(gen_leaves))),
                )
                # keep legacy fractional coverage too
                rec["coverage_at_100_frac"] = coverage_at_100(
                    g["leaves"], gen_leaves[:100] if len(gen_leaves) >= 100 else gen_leaves
                )
                rec["lineage_rf"] = float("nan")  # needs gen Newick export — TODO
                rec["meta"] = getattr(gt, "meta", {})
            except NotImplementedError as e:
                rec["status"] = f"not_implemented: {e}"
            except Exception as e:
                rec["status"] = f"error: {e}"
            rows.append(rec)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out = {
        "data": str(args.data),
        "methods": methods,
        "n_groups": len(groups),
        "rows": rows,
        "note": "Scaffold eval — Neutral SHM / AR stubs; CDR mask IMGT stub; root may be majority-vote.",
    }
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"n_rows": len(rows), "out": str(args.out)}, indent=2))


if __name__ == "__main__":
    main()
