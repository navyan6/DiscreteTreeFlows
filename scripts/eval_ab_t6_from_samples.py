#!/usr/bin/env python3
"""
Table 6 metrics from saved antibody_benchmark rollouts (Track C / Track B).

Computes per-method means over Rodriguez 82 families:
  Coverage@100 at absolute Hamming e∈{1,2,3,5}
  SHM load error
  CDR mutation recall (PCP codon coords → AA mask)
  Terminal diversity error (mean pairwise Hamming abs-diff)
  Lineage RF → NaN when topology forced or no gen Newick

Usage:
  python scripts/eval_ab_t6_from_samples.py \\
    --trees antibody_benchmark/data/processed/benchmark_trees.jsonl \\
    --samples-dir antibody_benchmark/results/samples \\
    --pcp ~/antibody_benchmark_raw/dasm/extracted/dasm-experiments-data/v3/rodriguez-....csv.gz \\
    --models thrifty,dasm_thrifty,cosine,treesbm,treesbm_ab \\
    --out benchmarks/results/tables/table6_ab_track_c.json
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from antibody_benchmark.codon import looks_like_nt, translate_nt  # noqa: E402
from scripts.ab_t6_metrics import (  # noqa: E402
    cdr_mask_from_pcp_codon_coords,
    cdr_mask_imgt_stub,
    cdr_mut_recall,
    coverage_at_k_abs,
    shm_load_error,
    terminal_diversity_error,
)


# Topology forced → RF not a free-generation metric
FORCED_TOPOLOGY = {
    "thrifty",
    "dasm_thrifty",
    "neutral_shm",
    "plm_prior",
    "treesbm",
    "treesbm_ab",
    "treesbm_ab_oas",
    "treesbm_ab_oas_v2",
    "treesbm_pathogen_ckpt",
    "identity_null",
    "poisson_null",
}


def _to_aa(seq: str) -> str:
    return translate_nt(seq) if looks_like_nt(seq) else seq


def _load_trees(path: Path) -> list[dict]:
    text = path.read_text().strip()
    if path.suffix == ".jsonl" or (text and text.startswith("{") and "\n{" in text + "\n"):
        raw = [json.loads(line) for line in text.splitlines() if line.strip()]
    else:
        raw = json.loads(text)
        if isinstance(raw, dict) and "trees" in raw:
            raw = raw["trees"]
    return raw


def _load_pcp_cdr_by_family(pcp_path: Path | None) -> dict[str, dict[str, tuple[int, int]]]:
    """family_id → {cdr1:(start,end), ...} from Rodriguez PCP (NT codon indices)."""
    out: dict[str, dict[str, tuple[int, int]]] = {}
    if pcp_path is None or not pcp_path.is_file():
        return out
    opener = gzip.open if str(pcp_path).endswith(".gz") else open
    with opener(pcp_path, "rt") as f:
        reader = csv.DictReader(f)
        for row in reader:
            fid = f"{row['sample_id']}::{row['family']}"
            if fid in out:
                continue
            coords = {}
            for k in (1, 2, 3):
                s = row.get(f"cdr{k}_codon_start_heavy")
                e = row.get(f"cdr{k}_codon_end_heavy")
                if s is None or e is None or s == "" or e == "":
                    continue
                try:
                    coords[f"cdr{k}"] = (int(float(s)), int(float(e)))
                except ValueError:
                    continue
            if coords:
                out[fid] = coords
    return out


def _leaf_aa(tree: dict) -> tuple[str, list[str]]:
    root_id = tree["root_id"]
    aa = tree.get("true_aa_sequences") or {}
    nt = tree.get("true_sequences") or {}
    root = aa.get(root_id) or tree.get("root_aa_sequence") or _to_aa(nt.get(root_id, tree.get("root_sequence", "")))
    leaves = []
    for lid in tree.get("leaf_ids") or [n for n, v in (tree.get("is_leaf") or {}).items() if v]:
        if lid == root_id:
            continue
        seq = aa.get(lid) or _to_aa(nt.get(lid, ""))
        if seq:
            leaves.append(seq)
    return root, leaves


def _gen_leaf_aa(
    tree: dict,
    seqs: dict[str, str],
    *,
    gen_root_id: str | None = None,
    free_topology: bool = False,
) -> list[str]:
    root_id = gen_root_id or tree["root_id"]
    leaf_ids = tree.get("leaf_ids") or [n for n, v in (tree.get("is_leaf") or {}).items() if v]
    out = []
    if not free_topology:
        for lid in leaf_ids:
            if lid == root_id:
                continue
            if lid not in seqs:
                continue
            out.append(_to_aa(seqs[lid]))
    # Free-topo / CoSiNE / mismatched ids — all non-root sequences as leaves
    if not out:
        for nid, s in seqs.items():
            if nid == root_id:
                continue
            # Skip obvious internal markers if any; keep all other nodes as pool
            out.append(_to_aa(s))
    return out


def _load_generated(samples_dir: Path, model: str) -> dict[str, list[dict]]:
    """family_id → list of {sequences, root_id?} payloads (sequences required)."""
    by_fam: dict[str, list[dict]] = defaultdict(list)
    model_dir = samples_dir / model
    if not model_dir.is_dir():
        return {}
    for path in sorted(model_dir.glob("*/rollout_*.json")):
        payload = json.loads(path.read_text())
        by_fam[payload["family_id"]].append(
            {
                "sequences": payload["sequences"],
                "root_id": payload.get("root_id"),
                "free_topology": bool(payload.get("free_topology", False)),
            }
        )
    return by_fam


def _mean(xs: list[float]) -> float:
    ys = [x for x in xs if x == x and not math.isinf(x)]
    return sum(ys) / len(ys) if ys else float("nan")


def _rf_note(model: str) -> tuple[float, str]:
    if model in FORCED_TOPOLOGY or model.startswith("treesbm"):
        return float("nan"), "topology_forced_RF_N/A"
    return float("nan"), "no_gen_newick_exported_RF_N/A"


def evaluate_model(
    trees: list[dict],
    gens: dict[str, list[dict[str, str]]],
    model: str,
    cdr_by_fam: dict[str, dict[str, tuple[int, int]]],
    *,
    k: int = 100,
    eps: list[int] | None = None,
    seed: int = 0,
) -> dict:
    eps = eps or [1, 2, 3, 5]
    rf_val, rf_note = _rf_note(model)
    per_family = []
    cov_lists = {e: [] for e in eps}
    shm_list, cdr_list, div_list = [], [], []

    for tree in trees:
        fid = tree["family_id"]
        root, gt_leaves = _leaf_aa(tree)
        rollouts = gens.get(fid, [])
        if not root or not gt_leaves or not rollouts:
            continue
        gen_leaves: list[str] = []
        for item in rollouts:
            if isinstance(item, dict) and "sequences" in item:
                seqs = item["sequences"]
                gen_leaves.extend(
                    _gen_leaf_aa(
                        tree,
                        seqs,
                        gen_root_id=item.get("root_id"),
                        free_topology=bool(item.get("free_topology", False)),
                    )
                )
            else:
                # legacy: bare sequences dict
                gen_leaves.extend(_gen_leaf_aa(tree, item))  # type: ignore[arg-type]
        if not gen_leaves:
            continue

        coords = cdr_by_fam.get(fid)
        mask = cdr_mask_from_pcp_codon_coords(len(root), coords) if coords else cdr_mask_imgt_stub(len(root))
        mask_source = "pcp" if coords else "imgt_stub"

        fam = {
            "family_id": fid,
            "n_gt_leaves": len(gt_leaves),
            "n_gen_leaves": len(gen_leaves),
            "n_rollouts": len(rollouts),
            "pool_k": min(k, len(gen_leaves)),
            "cdr_mask_source": mask_source,
            "cdr_mask_n_sites": int(sum(mask)),
        }
        for e in eps:
            c = coverage_at_k_abs(gt_leaves, gen_leaves, e=e, k=k, seed=seed + hash(fid) % 10_000)
            fam[f"coverage_at_{k}_e{e}"] = c
            cov_lists[e].append(c)
        shm = shm_load_error(root, gen_leaves, gt_leaves)
        cdr = cdr_mut_recall(root, gen_leaves, gt_leaves, mask)
        div = terminal_diversity_error(gen_leaves, gt_leaves)
        fam["shm_load_error"] = shm
        fam["cdr_mut_recall"] = cdr
        fam["terminal_diversity_error"] = div
        fam["lineage_rf"] = rf_val
        per_family.append(fam)
        shm_list.append(shm)
        cdr_list.append(cdr)
        div_list.append(div)

    summary = {
        "model": model,
        "n_families": len(per_family),
        "K": k,
        "lineage_rf": rf_val,
        "lineage_rf_note": rf_note,
        "shm_load_error": _mean(shm_list),
        "cdr_mut_recall": _mean(cdr_list),
        "terminal_diversity_error": _mean(div_list),
    }
    for e in eps:
        summary[f"coverage_at_{k}_e{e}"] = _mean(cov_lists[e])
    return {"summary": summary, "per_family": per_family}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", type=Path, default=ROOT / "antibody_benchmark/data/processed/benchmark_trees.jsonl")
    ap.add_argument("--samples-dir", type=Path, default=ROOT / "antibody_benchmark/results/samples")
    ap.add_argument(
        "--pcp",
        type=Path,
        default=None,
        help="Rodriguez PCP csv.gz for CDR codon coords (optional)",
    )
    ap.add_argument(
        "--models",
        default="neutral_shm,plm_prior,ar_tree_edit,thrifty,dasm_thrifty,cosine,treesbm,treesbm_ab,treesbm_ab_oas",
    )
    ap.add_argument("--K", type=int, default=100)
    ap.add_argument("--eps", default="1,2,3,5")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--md-out", type=Path, default=None)
    args = ap.parse_args()

    eps = [int(x) for x in args.eps.split(",") if x.strip()]
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    trees = _load_trees(args.trees)
    pcp = args.pcp
    if pcp is None:
        cand = Path.home() / (
            "antibody_benchmark_raw/dasm/extracted/dasm-experiments-data/v3/"
            "rodriguez-airr-seq-race-prod-NoWinCheck_igh_pcp_2024-11-12_MASKED_NI_noN_no-naive.csv.gz"
        )
        if cand.is_file():
            pcp = cand
        else:
            alt = ROOT / (
                "antibody_benchmark/data/raw/dasm/extracted/dasm-experiments-data/v3/"
                "rodriguez-airr-seq-race-prod-NoWinCheck_igh_pcp_2024-11-12_MASKED_NI_noN_no-naive.csv.gz"
            )
            pcp = alt if alt.is_file() else None
    cdr_by_fam = _load_pcp_cdr_by_family(pcp)
    print(f"trees={len(trees)} cdr_families={len(cdr_by_fam)} pcp={pcp}", flush=True)

    results = {}
    rows = []
    for model in models:
        gens = _load_generated(args.samples_dir, model)
        print(f"model={model} families_with_samples={len(gens)}", flush=True)
        if not gens:
            results[model] = {"summary": {"model": model, "n_families": 0, "status": "no_samples"}, "per_family": []}
            continue
        block = evaluate_model(trees, gens, model, cdr_by_fam, k=args.K, eps=eps, seed=args.seed)
        results[model] = block
        rows.append(block["summary"])
        s = block["summary"]
        print(
            f"  n={s['n_families']} cov@e2={s.get('coverage_at_100_e2')} "
            f"shm={s.get('shm_load_error')} cdr={s.get('cdr_mut_recall')} "
            f"tdiv={s.get('terminal_diversity_error')} rf={s.get('lineage_rf_note')}",
            flush=True,
        )

    payload = {
        "trees": str(args.trees),
        "samples_dir": str(args.samples_dir),
        "pcp": str(pcp) if pcp else None,
        "n_trees": len(trees),
        "n_cdr_masks": len(cdr_by_fam),
        "K": args.K,
        "eps": eps,
        "train_test_notes": {
            "neutral_shm": "JC69 NT CTMC on observed topo+BL (non-Thrifty Neutral SHM) → Rodriguez 82",
            "plm_prior": "ESM-2 8M mutation prior on observed topo+BL → Rodriguez 82",
            "ar_tree_edit": "ARTreeFormer N16 pool pruned to n_leaves + JTT (adapted; not Ab-trained) → Rodriguez 82",
            "thrifty": "ThriftyHumV0.2-59 (context SHM rates) → Rodriguez 82 — NOT paper Neutral SHM row",
            "dasm_thrifty": "DASMHumV1.0-4M + Thrifty → Rodriguez 82",
            "cosine": "CoSiNE dasm ckpt, unguided Gillespie on observed edges → Rodriguez 82",
            "treesbm": "pathogen checkpoints/best.pt (job 7062213) → Rodriguez 82 (OOD; paper Track C)",
            "treesbm_ab": "ab_dasm_v1 (Track B; Tang+VanWinkle export, job 7539619) → Rodriguez 82",
            "treesbm_ab_oas": "ab_oas_1m_v1 (Track A OAS 1M, CDR mask + entropy, λ_br=0) → Rodriguez 82",
            "treesbm_ab_oas_v2": "ab_oas_1m_v2_shm (Track A OAS 1M + SHM Q0 site-rate prior CDR∪AID) → Rodriguez 82",
        },
        "models": results,
        "table_rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n")

    md_path = args.md_out or args.out.with_suffix(".md")
    cols = [
        "model",
        "n_families",
        *[f"coverage_at_{args.K}_e{e}" for e in eps],
        "shm_load_error",
        "cdr_mut_recall",
        "terminal_diversity_error",
        "lineage_rf",
        "lineage_rf_note",
    ]

    def fmt(v):
        if v is None:
            return "—"
        if isinstance(v, float):
            if v != v:
                return "NaN"
            return f"{v:.4f}"
        return str(v)

    lines = [
        "# Table 6 — Antibody Track C / B (Rodriguez 82)",
        "",
        f"Source: `{args.out}`",
        f"CDR masks: PCP codon coords ({len(cdr_by_fam)} families)" if cdr_by_fam else "CDR masks: IMGT stub",
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
            "- **Coverage@100**: absolute Hamming ε; pool = up to K=100 gen leaf AAs across N=20 rollouts.",
            "- **Terminal diversity error**: |mean pairwise Hamming(gen) − mean pairwise Hamming(gt)| "
            "(scalar companion to leaf-div W1 in primary rollout table).",
            "- **Lineage RF**: NaN — TreeSBM/Thrifty/DASM force observed topology; CoSiNE samples lack Newick.",
            "- **treesbm** = pathogen ckpt; **treesbm_ab** = Track B `ab_dasm_v1`; **treesbm_ab_oas** = Track A `ab_oas_1m_v1`; **treesbm_ab_oas_v2** = Track A `ab_oas_1m_v2_shm` (SHM Q0 site-rate prior).",
            "- OAS v1 samples: `results/samples/treesbm_ab_oas`. v2: `results/samples/treesbm_ab_oas_v2` (do not overwrite `treesbm` / `treesbm_ab` / `treesbm_ab_oas`).",
            "",
        ]
    )
    md_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.out}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
