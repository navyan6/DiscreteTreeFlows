#!/usr/bin/env python3
"""
Table 5 — held-out *clonal lineage* split (locked).

Split unit = clone_id (donor|v_family|cdr3), NOT V-gene family and NOT random leaves.
Also enforces donor leakage ban: if a donor appears in test, none of that donor's
clones go to train (optional strict mode; default soft = clone-level only).

Usage:
  python scripts/ab_family_holdout_split.py \\
    --trees-dir data/ab_t5_500k/trees_all \\
    --clones-index data/ab_t5_500k/clones/clones_index.tsv \\
    --out-base data/ab_clones \\
    --test-frac 0.15 --val-frac 0.10 \\
    --donor-disjoint
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees-dir", type=Path, required=True)
    ap.add_argument("--clones-index", type=Path, required=True)
    ap.add_argument("--out-base", type=Path, default=Path("data/ab_clones"))
    ap.add_argument("--test-frac", type=float, default=0.15)
    ap.add_argument("--val-frac", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--donor-disjoint", action="store_true", default=True)
    ap.add_argument("--no-donor-disjoint", action="store_false", dest="donor_disjoint")
    args = ap.parse_args()

    rows = list(csv.DictReader(args.clones_index.open(), delimiter="\t"))
    # Map built groups to clone rows by order used in ab_build_clone_trees
    # (group_XXX built from sorted by_clone FASTAs matching safe_id).
    group_dirs = sorted(
        [p for p in args.trees_dir.iterdir() if p.is_dir() and p.name.startswith("group_")]
    )
    if not group_dirs:
        raise SystemExit(f"No group_* dirs in {args.trees_dir}")

    # Prefer matching via summary or safe_id in clone fasta name recorded in trees summary
    # Fallback: zip by sorted index order with clones_index sorted by -n_leaves (build used glob sort).
    by_safe = {r["safe_id"]: r for r in rows}
    assigned = []
    for gd in group_dirs:
        # find which clone: look for members tsv reference in sibling clones dir — use nwk parent name only
        # Heuristic: read leaf fasta first header comment — instead match via copies: original clone fasta path not stored.
        # Use trees summary.json list if present.
        assigned.append({"group_dir": gd, "meta": None})

    summary_path = args.trees_dir / "summary.json"
    # Rebuild mapping from by_clone sorted FASTA order (same as build script)
    clones_by = args.clones_index.parent / "by_clone"
    fastas = sorted(clones_by.glob("*.fasta")) if clones_by.is_dir() else []
    mapping = []
    for i, gd in enumerate(group_dirs):
        meta = None
        if i < len(fastas):
            safe = fastas[i].stem  # safe_id.fasta
            meta = by_safe.get(safe)
            if meta is None:
                # stem may equal safe_id
                meta = by_safe.get(fastas[i].name.replace(".fasta", ""))
        mapping.append({"group_dir": gd, "meta": meta, "safe_id": fastas[i].stem if i < len(fastas) else gd.name})

    rng = random.Random(args.seed)
    if args.donor_disjoint:
        donors = sorted({(m["meta"] or {}).get("donor_id", f"unk_{i}") for i, m in enumerate(mapping)})
        rng.shuffle(donors)
        n_test = max(1, int(len(donors) * args.test_frac))
        n_val = max(1, int(len(donors) * args.val_frac))
        test_d = set(donors[:n_test])
        val_d = set(donors[n_test : n_test + n_val])
        split_of = {}
        for m in mapping:
            d = (m["meta"] or {}).get("donor_id", "unknown")
            if d in test_d:
                split_of[m["group_dir"].name] = "test"
            elif d in val_d:
                split_of[m["group_dir"].name] = "val"
            else:
                split_of[m["group_dir"].name] = "train"
    else:
        idxs = list(range(len(mapping)))
        rng.shuffle(idxs)
        n_test = max(1, int(len(idxs) * args.test_frac))
        n_val = max(1, int(len(idxs) * args.val_frac))
        split_of = {}
        for j, i in enumerate(idxs):
            name = mapping[i]["group_dir"].name
            if j < n_test:
                split_of[name] = "test"
            elif j < n_test + n_val:
                split_of[name] = "val"
            else:
                split_of[name] = "train"

    counts = {"train": 0, "val": 0, "test": 0}
    for split in counts:
        (args.out_base / split).mkdir(parents=True, exist_ok=True)

    manifest = []
    for m in mapping:
        gd = m["group_dir"]
        split = split_of[gd.name]
        dest = args.out_base / split
        # flatten key artifacts to dest with renumbered group ids later
        new_idx = counts[split]
        counts[split] += 1
        gname = f"group_{new_idx:03d}"
        for src in gd.glob("*"):
            # normalize names
            suffix = src.name
            for old in (gd.name,):
                suffix = suffix.replace(old, gname)
            shutil.copy2(src, dest / suffix)
        # also expect group_XXX.fasta / .nwk at top for TreeDataset
        fa = gd / f"{gd.name}.fasta"
        nw = gd / f"{gd.name}.nwk"
        if fa.is_file():
            shutil.copy2(fa, dest / f"{gname}.fasta")
        if nw.is_file():
            shutil.copy2(nw, dest / f"{gname}.nwk")
        root = gd / f"{gd.name}_root.fasta"
        if root.is_file():
            shutil.copy2(root, dest / f"{gname}_root.fasta")
        manifest.append(
            {
                "split": split,
                "group": gname,
                "src_group": gd.name,
                "safe_id": m.get("safe_id"),
                "clone_id": (m["meta"] or {}).get("clone_id"),
                "donor_id": (m["meta"] or {}).get("donor_id"),
                "n_leaves": (m["meta"] or {}).get("n_leaves"),
            }
        )

    man_path = args.out_base / "split_manifest.tsv"
    with man_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["split", "group", "src_group", "safe_id", "clone_id", "donor_id", "n_leaves"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(manifest)

    summary = {
        "out_base": str(args.out_base),
        "holdout": "clonal_lineage",
        "donor_disjoint": args.donor_disjoint,
        "counts": counts,
        "n_mapped": len(mapping),
        "seed": args.seed,
    }
    (args.out_base / "split_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
