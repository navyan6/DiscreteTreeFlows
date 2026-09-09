#!/usr/bin/env python3
"""
Compare candidate root translations so the right repair is chosen per dataset.

`nt_to_aa` routes on `_frameshifting_gaps`: a sequence with a non-codon gap run
goes down the gap-aware, alignment-coordinate path. That routing is only valid
if the alignment is codon-framed, which `audit_frame_offset.py` shows is true
for the reference-anchored COVID spike and false for the flu and HIV
alignments, where no column offset gives a clean frame. Routing three quarters
of HIV down a path its alignment cannot support would trade a truncation bug
for a garbage-protein bug.

Three candidates per root:

  legacy   strip gaps, translate, truncate at the first stop  (current, buggy)
  notrunc  strip gaps, translate, keep full length            (stops -> X)
  aligned  translate in alignment coordinates                 (the COVID repair)

Judged against the modal leaf length in the same tree. The winner is the one
that restores the root to leaf length without inventing stop codons.
"""

from __future__ import annotations

import argparse
import collections
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def is_root(name: str) -> bool:
    return name == "root" or name.endswith("NODE_0000000")


def is_leaf(name: str) -> bool:
    if "|" in name:
        return name.rsplit("|", 1)[1] == "leaf"
    return not name.startswith("NODE_") and name != "root"


def _tr(cds: str) -> str:
    cds = cds[:len(cds) - len(cds) % 3]
    out = []
    for i in range(0, len(cds), 3):
        c = cds[i:i + 3]
        try:
            out.append(str(Seq(c).translate()))
        except Exception:  # noqa: BLE001
            out.append("X")
    return "".join(out)


def legacy(aln: str) -> str:
    s = aln.replace("-", "")
    i = s.find("ATG")
    if i == -1:
        return ""
    aa = _tr(s[i:])
    k = aa.find("*")
    return aa[:k] if k != -1 else aa


def notrunc(aln: str) -> str:
    s = aln.replace("-", "")
    i = s.find("ATG")
    if i == -1:
        return ""
    aa = _tr(s[i:])
    return aa.rstrip("*").replace("*", "X")


def aligned(aln: str) -> str:
    ung = aln.replace("-", "")
    i = ung.find("ATG")
    if i == -1:
        return ""
    seen, start = 0, None
    for k, ch in enumerate(aln):
        if ch == "-":
            continue
        if seen == i:
            start = k
            break
        seen += 1
    if start is None:
        return ""
    out = []
    for j in range(start, len(aln) - 2, 3):
        c = aln[j:j + 3]
        g = c.count("-")
        if g == 3:
            continue
        if g:
            out.append("X")
        else:
            try:
                out.append(str(Seq(c).translate()))
            except Exception:  # noqa: BLE001
                out.append("X")
    return "".join(out).rstrip("*").replace("*", "X")


def run(d: Path, max_trees: int) -> None:
    files = sorted(d.glob("group_*_anc_nt.fasta"))[:max_trees]
    name = "/".join(str(d).split("/")[-2:])
    if not files:
        print(f"\n=== {name} === (none)")
        return

    stats = {k: {"ratio": [], "stops": 0, "ok": 0} for k in ("legacy", "notrunc", "aligned")}
    n = 0
    for f in files:
        recs = list(SeqIO.parse(f, "fasta"))
        root = next((r for r in recs if is_root(r.id)), None)
        leaves = [len(str(r.seq).replace("-", "")) for r in recs if is_leaf(r.id)]
        if root is None or not leaves:
            continue
        modal_nt = collections.Counter(leaves).most_common(1)[0][0]
        target = modal_nt // 3
        s = str(root.seq).upper()
        n += 1
        for key, fn in (("legacy", legacy), ("notrunc", notrunc), ("aligned", aligned)):
            aa = fn(s)
            r = len(aa) / target if target else 0.0
            stats[key]["ratio"].append(r)
            stats[key]["stops"] += aa.count("X")
            if 0.95 <= r <= 1.05:
                stats[key]["ok"] += 1

    print(f"\n=== {name} ===  {n} trees")
    print(f"  {'candidate':<10}{'median len/leaf':>17}{'at leaf length':>16}{'X per root':>12}")
    for key in ("legacy", "notrunc", "aligned"):
        rs = sorted(stats[key]["ratio"])
        med = rs[len(rs) // 2] if rs else 0
        print(f"  {key:<10}{med:>17.3f}{stats[key]['ok']:>10}/{n:<5}"
              f"{stats[key]['stops'] / max(n, 1):>12.1f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--max-trees", type=int, default=40)
    a = ap.parse_args()
    for d in a.dirs:
        run(d, a.max_trees)


if __name__ == "__main__":
    main()
