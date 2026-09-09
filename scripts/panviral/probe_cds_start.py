#!/usr/bin/env python3
"""
Check the assumption that the CDS begins at the first ATG.

`nt_to_aa` locates the start codon with `seq.find('ATG')`. That is safe only if
the alignment window begins at the true start codon, because the first ATG in an
arbitrary window is frequently an internal methionine. If its index is not a
multiple of three the whole translation is read out of frame and dies at the
first in-frame stop, which is how Bundibugyo ended up with 17-residue proteins
against a 900-codon alignment that has no gaps at all.

Prints, per dataset, where the first ATG sits and what the protein looks like
when translated from column 0 instead.
"""

from __future__ import annotations

import argparse
import collections
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def is_leaf(name: str) -> bool:
    if "|" in name:
        return name.rsplit("|", 1)[1] == "leaf"
    return not name.startswith("NODE_") and name != "root"


def translate_from(aln: str, start: int) -> str:
    out = []
    for i in range(start, len(aln) - 2, 3):
        c = aln[i:i + 3]
        if c == "---":
            continue
        if "-" in c:
            out.append("X")
            continue
        try:
            out.append(str(Seq(c).translate()))
        except Exception:  # noqa: BLE001
            out.append("X")
    return "".join(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dirs", nargs="+", type=Path)
    a = ap.parse_args()

    hdr = (f"{'dataset':<34}{'starts ATG':>11}{'1st ATG':>9}{'%3':>4}"
           f"{'len@ATG':>9}{'len@col0':>10}{'stops@0':>9}")
    print(hdr)
    print("-" * len(hdr))
    for d in a.dirs:
        files = sorted(d.glob("group_*_anc_nt.fasta"))[:5]
        name = "/".join(str(d).split("/")[-2:])
        if not files:
            print(f"{name[-33:]:<34} (none)")
            continue
        atg_idx: collections.Counter = collections.Counter()
        starts = n = 0
        len_atg, len_col0, stops0 = [], [], []
        for f in files:
            for r in SeqIO.parse(f, "fasta"):
                if not is_leaf(r.id):
                    continue
                s = str(r.seq).upper()
                n += 1
                if s.startswith("ATG"):
                    starts += 1
                i = s.replace("-", "").find("ATG")
                atg_idx[i] += 1
                p_atg = translate_from(s.replace("-", ""), max(i, 0))
                k = p_atg.find("*")
                len_atg.append(len(p_atg[:k] if k != -1 else p_atg))
                p0 = translate_from(s, 0)
                len_col0.append(len(p0.rstrip("*")))
                stops0.append(p0.rstrip("*").count("*"))
                if n >= 40:
                    break
        if not n:
            print(f"{name[-33:]:<34} (no leaves)")
            continue
        mi = atg_idx.most_common(1)[0][0]
        med = lambda v: sorted(v)[len(v) // 2]  # noqa: E731
        print(f"{name[-33:]:<34}{f'{starts}/{n}':>11}{mi:>9}{mi % 3:>4}"
              f"{med(len_atg):>9}{med(len_col0):>10}{med(stops0):>9}")


if __name__ == "__main__":
    main()
