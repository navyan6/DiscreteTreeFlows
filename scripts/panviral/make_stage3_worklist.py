#!/usr/bin/env python3
"""
Build the stage-3 worklist: every virus split that has input groups and still
needs trees.

One line per (virus, split): slug<TAB>split<TAB>n_groups<TAB>data_dir
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path("data/panviral"))
    ap.add_argument("--out", type=Path, default=Path("data/panviral/stage3_worklist.tsv"))
    ap.add_argument("--splits", nargs="+", default=["train", "test"])
    ap.add_argument("--redo", action="store_true",
                    help="include splits that already have rooted trees")
    args = ap.parse_args()

    rows = []
    for manifest in sorted(args.root.glob("*/manifest.json")):
        slug = manifest.parent.name
        try:
            info = json.loads(manifest.read_text())
        except Exception:  # noqa: BLE001
            continue
        for split in args.splits:
            d = manifest.parent / split
            if not d.is_dir():
                continue
            groups = sorted(d.glob(f"{slug}_group_*.fasta"))
            if not groups:
                continue
            rooted = list(d.glob("group_*_rooted.nwk"))
            if rooted and not args.redo and len(rooted) >= len(groups):
                continue
            rows.append((slug, split, len(groups), str(d)))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {len(rows)} split(s) -> {args.out}")


if __name__ == "__main__":
    main()
