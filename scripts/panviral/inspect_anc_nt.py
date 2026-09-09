#!/usr/bin/env python3
"""
Characterise ancestral nucleotide files before retranslating anything.

The COVID repair is not automatically the right repair elsewhere. These are
different genes, aligned separately, and the gap-aware translator only helps
where the old one actually frameshifted -- that is, where a sequence carries a
gap run whose length is not a multiple of three. A dataset whose gaps are all
codon-sized was translated correctly the first time and must come out
byte-identical; running a "fix" there would be churn at best.

Reports per dataset: alignment width and whether it is codon-commensurate, how
many sequences carry gaps at all, how many carry frameshifting gaps, and how
many roots are short relative to their leaves.
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from Bio import SeqIO  # noqa: E402

try:
    from treeencoder.seq_utils import _frameshifting_gaps  # noqa: E402
except Exception:  # noqa: BLE001 - fall back to a local copy of the rule
    import re

    def _frameshifting_gaps(seq: str) -> bool:
        return any(len(m.group(0)) % 3 for m in re.finditer(r"-+", seq))


def is_root(name: str) -> bool:
    return name == "root" or name.endswith("NODE_0000000")


def is_leaf(name: str) -> bool:
    if "|" in name:
        return name.rsplit("|", 1)[1] == "leaf"
    return not name.startswith("NODE_") and name != "root"


def inspect(d: Path) -> dict | None:
    files = sorted(d.glob("group_*_anc_nt.fasta"))
    if not files:
        return None

    widths: collections.Counter = collections.Counter()
    n_seq = n_gapped = n_frameshift = 0
    n_fs_root = n_fs_leaf = 0
    short_roots = 0
    trees_with_fs = 0

    for f in files:
        recs = list(SeqIO.parse(f, "fasta"))
        if not recs:
            continue
        widths[len(recs[0].seq)] += 1
        tree_fs = False
        leaf_lens = []
        root_len = None
        for r in recs:
            s = str(r.seq).upper()
            n_seq += 1
            if "-" in s:
                n_gapped += 1
                if _frameshifting_gaps(s):
                    n_frameshift += 1
                    tree_fs = True
                    if is_root(r.id):
                        n_fs_root += 1
                    elif is_leaf(r.id):
                        n_fs_leaf += 1
            ungapped = len(s.replace("-", ""))
            if is_root(r.id):
                root_len = ungapped
            elif is_leaf(r.id):
                leaf_lens.append(ungapped)
        if tree_fs:
            trees_with_fs += 1
        if root_len is not None and leaf_lens:
            modal_leaf = collections.Counter(leaf_lens).most_common(1)[0][0]
            if root_len < 0.9 * modal_leaf:
                short_roots += 1

    width, _ = widths.most_common(1)[0]
    return {
        "dir": str(d), "trees": len(files), "width": width,
        "width_mod3": width % 3, "seqs": n_seq,
        "gapped": n_gapped, "frameshifting": n_frameshift,
        "fs_roots": n_fs_root, "fs_leaves": n_fs_leaf,
        "trees_with_fs": trees_with_fs, "short_roots": short_roots,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dirs", nargs="+", type=Path)
    args = ap.parse_args()

    hdr = (f"{'dataset':<34}{'trees':>6}{'width':>7}{'%3':>4}{'seqs':>8}"
           f"{'gapped':>8}{'frameshift':>11}{'fs-trees':>9}{'short roots':>12}")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for d in args.dirs:
        r = inspect(d)
        if r is None:
            print(f"{str(d)[-33:]:<34}  (no group_*_anc_nt.fasta)")
            continue
        rows.append(r)
        name = r["dir"].split("/")[-2] + "/" + r["dir"].split("/")[-1]
        print(f"{name[-33:]:<34}{r['trees']:>6}{r['width']:>7}{r['width_mod3']:>4}"
              f"{r['seqs']:>8}{r['gapped']:>8}{r['frameshifting']:>11}"
              f"{r['trees_with_fs']:>9}{r['short_roots']:>12}")

    print("\nverdict")
    for r in rows:
        name = r["dir"].split("/")[-2] + "/" + r["dir"].split("/")[-1]
        if r["frameshifting"] == 0:
            v = "NO-OP -- no frameshifting gaps; retranslation must change nothing"
        else:
            pct = 100 * r["frameshifting"] / max(r["seqs"], 1)
            v = (f"REPAIR -- {r['frameshifting']} seqs ({pct:.1f}%) frameshifted "
                 f"across {r['trees_with_fs']} trees "
                 f"[{r['fs_roots']} roots, {r['fs_leaves']} leaves]")
        print(f"  {name:<34} {v}")
        if r["width_mod3"]:
            print(f"  {'':34} NOTE alignment width {r['width']} is not a multiple "
                  f"of 3 -- check the CDS window before trusting any translation")


if __name__ == "__main__":
    main()
