#!/usr/bin/env python3
"""
Separate stale protein files from malformed groups.

The retranslation guard refuses any change it cannot explain, and in
h1n1_epidemic/train it refused 106 groups. Two very different things were
hiding behind that count:

  stale     the stored _anc_aa file was not derived from the stored _anc_nt
            file. Re-running the ancestral reconstruction rewrote the
            nucleotides without rewriting the protein, so the two disagree at a
            handful of residues. Nucleotides are upstream, so retranslating is
            the fix and the disagreement is the evidence it was needed.

  malformed the nucleotide file itself is not a gene -- group_052 is 468 nt for
            an HA that should be ~1734. No translation rescues that; the group
            should be dropped rather than repaired.

Telling them apart matters because the first wants writing and the second wants
excluding, and the refusal count alone does not distinguish them.
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from Bio import SeqIO  # noqa: E402
from treeencoder.seq_utils import nt_to_aa  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--mode", default="ungapped")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="nt width below this fraction of modal = malformed")
    a = ap.parse_args()

    for d in a.dirs:
        files = sorted(d.glob("group_*_anc_nt.fasta"))
        name = "/".join(str(d).split("/")[-2:])
        if not files:
            print(f"\n=== {name} === (none)")
            continue

        widths = {}
        for f in files:
            r = next(SeqIO.parse(f, "fasta"), None)
            if r is not None:
                widths[f] = len(r.seq)
        modal = collections.Counter(widths.values()).most_common(1)[0][0]

        malformed, stale, clean = [], [], []
        for f in files:
            w = widths.get(f, 0)
            if w < a.min_frac * modal:
                malformed.append((f.name, w))
                continue
            aa_path = f.with_name(f.name.replace("_anc_nt", "_anc_aa"))
            if not aa_path.exists():
                continue
            old = {r.id: str(r.seq) for r in SeqIO.parse(aa_path, "fasta")}
            n_bad = n_tot = 0
            for r in SeqIO.parse(f, "fasta"):
                o = old.get(r.id)
                if o is None:
                    continue
                n_tot += 1
                try:
                    n = nt_to_aa(str(r.seq).upper(), mode=a.mode)
                except Exception:  # noqa: BLE001
                    n_bad += 1
                    continue
                if not (n.startswith(o) or o.startswith(n)):
                    n_bad += 1
            (stale if n_bad else clean).append((f.name, n_bad, n_tot))

        print(f"\n=== {name} ===  {len(files)} groups, modal nt width {modal}")
        print(f"  clean      {len(clean):>4}  (old protein is a prefix of new)")
        print(f"  stale      {len(stale):>4}  (protein file out of sync with nucleotides)")
        print(f"  malformed  {len(malformed):>4}  (nucleotide file far below modal width)")
        for nm, w in malformed[:6]:
            print(f"      {nm}  {w} nt  ({w / modal:.0%} of modal)")
        if stale:
            tot_bad = sum(b for _, b, _ in stale)
            tot_seq = sum(t for _, _, t in stale)
            print(f"      {tot_bad} of {tot_seq} sequences disagree "
                  f"({tot_bad / max(tot_seq, 1):.1%}) across stale groups")


if __name__ == "__main__":
    main()
