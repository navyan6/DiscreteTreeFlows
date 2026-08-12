#!/usr/bin/env python3
"""Build and freeze held-out AntibodyTree artifacts from DASM PCP tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antibody_benchmark.rollout.reconstruct_trees import (
    format_attrition_summary,
    reconstruct_trees_from_pcp_df,
    summarize_trees,
)


def _resolve_raw_dasm(cfg: dict) -> Path:
    env = os.environ.get("AB_DASM_DIR") or os.environ.get("ANTIBODY_DASM_DIR")
    if env:
        return Path(env).expanduser().resolve()
    raw = cfg["paths"]["raw_dasm"]
    p = Path(raw).expanduser()
    if not p.is_absolute():
        p = ROOT / p
    return p.resolve()


def _file_fingerprint(path: Path) -> dict:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "md5": h.hexdigest(),
    }


def _tree_to_freeze_dict(t) -> dict:
    """Frozen family record. true_sequences are eval-only (not for rollout)."""
    return {
        "family_id": t.family_id,
        "donor_id": t.donor_id,
        "root_id": t.root_id,
        "root_sequence": t.true_sequences[t.root_id],
        "root_aa_sequence": t.true_aa_sequences.get(t.root_id, ""),
        "node_ids": t.node_ids,
        "leaf_ids": t.leaves(),
        "edges": [list(e) for e in t.edges],
        "branch_lengths": {f"{a}->{b}": float(v) for (a, b), v in t.branch_lengths.items()},
        "true_sequences": t.true_sequences,  # eval-only
        "true_aa_sequences": t.true_aa_sequences,  # eval-only
        "is_leaf": t.is_leaf,
        "metadata": {k: v for k, v in t.metadata.items() if k != "light_seqs"},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="antibody_benchmark/configs/default.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text())
    paths = cfg["paths"]
    ds = cfg["dataset"]

    raw = _resolve_raw_dasm(cfg)
    pcp = raw / ds["primary_pcp"]
    if not pcp.is_file():
        raise SystemExit(
            f"BLOCKER: PCP not found: {pcp}. "
            "Set AB_DASM_DIR to Betty/local extracted dasm-experiments-data "
            "(e.g. ~/antibody_benchmark_raw/dasm/extracted/dasm-experiments-data)."
        )

    fp = _file_fingerprint(pcp)
    df = pd.read_csv(pcp)
    trees, attrition_df, reason_counts, stage_counts = reconstruct_trees_from_pcp_df(
        df,
        chain=ds.get("chain", "heavy"),
        require_single_root=bool(ds.get("require_single_root", True)),
        min_leaves=int(ds.get("min_leaves", 4)),
        min_depth=int(ds.get("min_depth", 2)),
        translate_aa=True,
        max_families=ds.get("max_families"),
        return_attrition=True,
    )
    summary = summarize_trees(trees)

    out_proc = ROOT / paths["processed"]
    out_proc.mkdir(parents=True, exist_ok=True)
    if "splits" in paths:
        out_split = ROOT / paths["splits"]
        out_split.mkdir(parents=True, exist_ok=True)
    else:
        out_split = out_proc

    # Frozen primary artifacts
    jsonl_path = out_proc / "benchmark_trees.jsonl"
    with jsonl_path.open("w") as f:
        for t in trees:
            f.write(json.dumps(_tree_to_freeze_dict(t), sort_keys=True) + "\n")

    fam_path = out_proc / "benchmark_family_ids.txt"
    fam_path.write_text("\n".join(t.family_id for t in trees) + ("\n" if trees else ""))

    audit_csv = out_proc / "data_audit.csv"
    attrition_df.to_csv(audit_csv, index=False)

    filter_md = out_proc / "FILTER_SUMMARY.md"
    filter_md.write_text(
        format_attrition_summary(
            stage_counts,
            reason_counts,
            source_pcp=str(pcp),
            n_final=len(trees),
        )
        + "\n## Source fingerprint\n\n"
        + f"- path: `{fp['path']}`\n"
        + f"- size_bytes: {fp['size_bytes']}\n"
        + f"- md5: `{fp['md5']}`\n"
        + f"- AB_DASM_DIR: `{os.environ.get('AB_DASM_DIR', '')}`\n"
        + f"- root_source_counts: `{summary.get('root_source_counts', {})}`\n"
    )

    # Back-compat copies used by older scripts
    legacy = out_proc / "heldout_trees.json"
    legacy.write_text(json.dumps([_tree_to_freeze_dict(t) for t in trees], indent=2))
    (out_split / "heldout_families.txt").write_text(fam_path.read_text())
    (out_proc / "heldout_summary.json").write_text(
        json.dumps({"n_trees": len(trees), **summary, "source": fp}, indent=2)
    )

    print(
        json.dumps(
            {
                "n_trees": len(trees),
                "stage_counts": stage_counts,
                "reason_counts": reason_counts,
                "source": fp,
                "artifacts": {
                    "benchmark_trees.jsonl": str(jsonl_path),
                    "benchmark_family_ids.txt": str(fam_path),
                    "data_audit.csv": str(audit_csv),
                    "FILTER_SUMMARY.md": str(filter_md),
                },
                **summary,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
