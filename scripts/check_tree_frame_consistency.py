#!/usr/bin/env python3
"""
Is column-wise leaf-vs-root comparison safe inside a single tree?

eval_evescape_enrichment.py calls mutations by walking leaf and root together
column by column. That is only meaningful if the two sequences are in the same
frame. 74% of COVID leaves are not 1273 aa (NUMBERS_THAT_CHANGE.md), so the
question is whether a tree's leaves at least share their *root's* frame.

Equal length is necessary but not sufficient: two sequences can both be 1271 aa
via different deletions. So we also check Hamming distance for the equal-length
pairs. Same frame means a handful of substitutions; a frame shift smears the
whole downstream sequence and produces hundreds of mismatches.

  python scripts/check_tree_frame_consistency.py --data-dir data/covid/train
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

SHIFT_HAMMING = 50  # more mismatches than this against the root implies a frame shift


def iter_fasta(path: Path):
    name = None
    buf: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def is_leaf(name: str) -> bool:
    return not name.startswith("NODE_") and name not in ("root", "root|root")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=ROOT / "data" / "covid" / "train")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results" / "voc_threat_panel" / "tree_frame_consistency.json")
    args = ap.parse_args()

    trees = []
    tot_leaves = tot_same_len = tot_same_frame = 0
    ham_all: list[int] = []

    paths = sorted(args.data_dir.glob("group_*_anc_aa.fasta"))
    print(f"scanning {len(paths)} trees in {args.data_dir} …")
    for n, path in enumerate(paths, 1):
        recs = list(iter_fasta(path))
        internals = {k: v for k, v in recs if not is_leaf(k)}
        leaves = [v for k, v in recs if is_leaf(k)]
        if not leaves or not internals:
            continue
        root = internals.get("NODE_0000000") or next(iter(internals.values()))

        same_len = 0
        same_frame = 0
        hams: list[int] = []
        for s in leaves:
            if len(s) != len(root):
                continue
            same_len += 1
            h = sum(1 for a, b in zip(s, root) if a != b)
            hams.append(h)
            if h <= SHIFT_HAMMING:
                same_frame += 1

        tot_leaves += len(leaves)
        tot_same_len += same_len
        tot_same_frame += same_frame
        ham_all.extend(hams)
        trees.append(
            {
                "tree": path.name,
                "n_leaves": len(leaves),
                "root_len": len(root),
                "frac_same_len": round(same_len / len(leaves), 4),
                "frac_same_frame": round(same_frame / len(leaves), 4),
                "median_hamming_same_len": (
                    round(statistics.median(hams), 1) if hams else None
                ),
            }
        )
        if n % 100 == 0:
            print(f"  {n}/{len(paths)}")

    print(f"\ntrees: {len(trees)}   leaves: {tot_leaves}")
    print(f"  leaf length == root length      : {tot_same_len} ({tot_same_len / tot_leaves:.3f})")
    print(f"  ... and Hamming <= {SHIFT_HAMMING} vs root : {tot_same_frame} "
          f"({tot_same_frame / tot_leaves:.3f})")
    if ham_all:
        ham_all.sort()
        q = lambda f: ham_all[int(f * (len(ham_all) - 1))]  # noqa: E731
        print(f"  Hamming(leaf, root) on equal-length pairs: "
              f"median={q(0.5)} p90={q(0.9)} p99={q(0.99)} max={ham_all[-1]}")
        print(f"  equal-length pairs that still look frame-shifted (>{SHIFT_HAMMING}): "
              f"{sum(1 for h in ham_all if h > SHIFT_HAMMING)} / {len(ham_all)}")

    fr = Counter()
    for t in trees:
        fr[round(t["frac_same_frame"], 1)] += 1
    print("\n  per-tree frac of leaves sharing the root's frame:")
    for k in sorted(fr):
        print(f"    {k:>4.1f}: {fr[k]:>4} trees")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(
            {
                "data_dir": str(args.data_dir),
                "n_trees": len(trees),
                "n_leaves": tot_leaves,
                "frac_leaf_len_eq_root": round(tot_same_len / max(tot_leaves, 1), 4),
                "frac_leaf_frame_eq_root": round(tot_same_frame / max(tot_leaves, 1), 4),
                "trees": trees,
            },
            indent=2,
        )
        + "\n"
    )
    try:
        print(f"\nwrote {args.out.relative_to(ROOT)}")
    except ValueError:
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
