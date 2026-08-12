#!/usr/bin/env python3
"""Export DASM heavy PCP trees → TreeDataset layout for Ab TreeSBM training.

Uses Tang + VanWinkle-IGH train PCPs (excludes Rodriguez held-out). Writes:
  data/ab_dasm_trees/{train,val}/group_XXX_{rooted.nwk,anc_aa.fasta,bl.json}

Does NOT train on the 82 Rodriguez benchmark families.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antibody_benchmark.rollout.reconstruct_trees import reconstruct_trees_from_pcp_df

# Heavy-only DASM train sources (Rodriguez = held-out, never train).
DEFAULT_PCPS = [
    "v3/tang-deepshm-prod-NoWinCheck_igh_pcp_2024-10-29_MASKED_NI_ConsCys_no-naive_DXSMVALID.csv.gz",
    "v3/v3convert_vanwinkle-170-igh_pcp_2025-03-05_MASKED_NI_train_no-naive_DXSMVALID_ConsCys.csv.gz",
]

_SAFE = re.compile(r"[^A-Za-z0-9_.-]+")


def _sanitize(name: str) -> str:
    s = _SAFE.sub("_", str(name)).strip("_")
    return s or "NODE"


def _to_newick(tree) -> str:
    children: dict[str, list[str]] = defaultdict(list)
    for p, c in tree.edges:
        children[p].append(c)
    id_map = {n: _sanitize(n) for n in tree.node_ids}
    # Ensure uniqueness after sanitize
    seen: dict[str, int] = {}
    for n, s in list(id_map.items()):
        if s in seen:
            seen[s] += 1
            id_map[n] = f"{s}_{seen[s]}"
        else:
            seen[s] = 0

    def rec(nid: str, bl: float | None) -> str:
        kids = children.get(nid, [])
        label = id_map[nid]
        if not kids:
            return f"{label}:{float(bl):.8g}" if bl is not None else label
        parts = []
        for c in kids:
            parts.append(rec(c, float(tree.branch_lengths[(nid, c)])))
        inner = f"({','.join(parts)}){label}"
        if bl is None:
            return inner
        return f"{inner}:{float(bl):.8g}"

    return rec(tree.root_id, None) + ";"


def _write_group(out_dir: Path, g: int, tree, id_map: dict[str, str]) -> int:
    """Write one TreeDataset group. Returns AA seq length."""
    nwk = out_dir / f"group_{g:03d}_rooted.nwk"
    fa = out_dir / f"group_{g:03d}_anc_aa.fasta"
    blj = out_dir / f"group_{g:03d}_bl.json"

    # Rebuild newick with same sanitize map as fasta ids
    children: dict[str, list[str]] = defaultdict(list)
    for p, c in tree.edges:
        children[p].append(c)

    def rec(nid: str, bl: float | None) -> str:
        kids = children.get(nid, [])
        label = id_map[nid]
        if not kids:
            return f"{label}:{float(bl):.8g}" if bl is not None else label
        parts = [rec(c, float(tree.branch_lengths[(nid, c)])) for c in kids]
        inner = f"({','.join(parts)}){label}"
        return inner if bl is None else f"{inner}:{float(bl):.8g}"

    nwk.write_text(rec(tree.root_id, None) + ";\n")

    aa = tree.true_aa_sequences
    if not aa:
        raise RuntimeError(f"missing AA sequences for {tree.family_id}")
    lens = {len(s) for s in aa.values() if s}
    if len(lens) != 1:
        # pad/truncate to modal length of root
        L = len(aa[tree.root_id])
    else:
        L = next(iter(lens))

    lines = []
    for nid in tree.node_ids:
        seq = aa.get(nid, "")
        if not seq:
            continue
        seq = (seq + ("-" * L))[:L]
        lines.append(f">{id_map[nid]}\n{seq}")
    fa.write_text("\n".join(lines) + "\n")

    nodes = {id_map[n]: {"numdate": 0.0, "branch_length": 0.0} for n in tree.node_ids}
    for (p, c), bl in tree.branch_lengths.items():
        nodes[id_map[c]]["branch_length"] = float(bl)
    blj.write_text(
        json.dumps(
            {
                "generated_by": "antibody_benchmark.export_treesbm_train_data",
                "family_id": tree.family_id,
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n"
    )
    return L


def _id_map(tree) -> dict[str, str]:
    raw = {n: _sanitize(n) for n in tree.node_ids}
    seen: dict[str, int] = {}
    out: dict[str, str] = {}
    for n, s in raw.items():
        if s in seen:
            seen[s] += 1
            out[n] = f"{s}_{seen[s]}"
        else:
            seen[s] = 0
            out[n] = s
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--dasm-dir",
        default=None,
        help="Extracted dasm-experiments-data root (or set AB_DASM_DIR)",
    )
    ap.add_argument(
        "--out",
        default="data/ab_dasm_trees",
        help="Output base (train/val subdirs created)",
    )
    ap.add_argument("--min-leaves", type=int, default=4)
    ap.add_argument("--min-depth", type=int, default=2)
    ap.add_argument("--val-frac", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument(
        "--max-trees",
        type=int,
        default=800,
        help="Cap total trees (TreeDataset group_%%03d → max 999)",
    )
    ap.add_argument(
        "--pcps",
        nargs="*",
        default=DEFAULT_PCPS,
        help="PCP paths relative to dasm-dir",
    )
    args = ap.parse_args()

    import os

    dasm = Path(
        args.dasm_dir
        or os.environ.get("AB_DASM_DIR")
        or os.environ.get("ANTIBODY_DASM_DIR")
        or (ROOT / "antibody_benchmark/data/raw/dasm/extracted/dasm-experiments-data")
    ).expanduser().resolve()
    if not dasm.is_dir():
        raise SystemExit(f"BLOCKER: DASM dir not found: {dasm}")

    trees = []
    for rel in args.pcps:
        pcp = dasm / rel
        if not pcp.is_file():
            raise SystemExit(f"BLOCKER: missing PCP {pcp}")
        if "rodriguez" in pcp.name.lower():
            raise SystemExit(f"REFUSING held-out PCP: {pcp.name}")
        print(f"Loading {pcp.name} ...", flush=True)
        df = pd.read_csv(pcp)
        got, *_ = reconstruct_trees_from_pcp_df(
            df,
            chain="heavy",
            require_single_root=True,
            min_leaves=args.min_leaves,
            min_depth=args.min_depth,
            translate_aa=True,
            return_attrition=True,
        )
        print(f"  eligible: {len(got)}", flush=True)
        trees.extend(got)

    # Dedup by family_id
    by_id = {t.family_id: t for t in trees}
    trees = list(by_id.values())
    rng = __import__("random").Random(args.seed)
    rng.shuffle(trees)
    if args.max_trees and len(trees) > args.max_trees:
        trees = trees[: args.max_trees]
    n_val = max(1, int(round(len(trees) * args.val_frac))) if len(trees) > 1 else 0
    val_trees = trees[:n_val]
    train_trees = trees[n_val:]

    out_base = Path(args.out)
    if not out_base.is_absolute():
        out_base = ROOT / out_base
    meta = {
        "dasm_dir": str(dasm),
        "pcps": list(args.pcps),
        "n_train": len(train_trees),
        "n_val": len(val_trees),
        "min_leaves": args.min_leaves,
        "min_depth": args.min_depth,
        "seed": args.seed,
        "held_out_excluded": "rodriguez",
        "max_aa_len": 0,
    }

    max_L = 0
    for split, subset in (("train", train_trees), ("val", val_trees)):
        d = out_base / split
        d.mkdir(parents=True, exist_ok=True)
        # clear stale groups
        for pat in ("group_*_rooted.nwk", "group_*_anc_aa.fasta", "group_*_bl.json"):
            for old in d.glob(pat):
                old.unlink()
        for i, tree in enumerate(subset):
            imap = _id_map(tree)
            L = _write_group(d, i, tree, imap)
            max_L = max(max_L, L)
            if (i + 1) % 50 == 0 or i == 0:
                print(f"  wrote {split} {i}/{len(subset)} L={L}", flush=True)
        print(f"Wrote {len(subset)} trees → {d}", flush=True)

    meta["max_aa_len"] = max_L
    meta_path = out_base / "export_meta.json"
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    print(f"Recommend --max-seq-len >= {max_L} (use 200 if max_aa_len<=200)")


if __name__ == "__main__":
    main()
