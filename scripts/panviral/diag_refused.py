#!/usr/bin/env python3
"""Show why a sequence's retranslation was refused: old vs new, side by side."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from Bio import SeqIO  # noqa: E402
from treeencoder.seq_utils import nt_to_aa  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--group", type=Path, required=True, help="path to *_anc_nt.fasta")
    ap.add_argument("--name", default="NODE_0000000")
    ap.add_argument("--mode", default="ungapped")
    a = ap.parse_args()

    nt = {r.id: str(r.seq).upper() for r in SeqIO.parse(a.group, "fasta")}
    aa_path = a.group.with_name(a.group.name.replace("_anc_nt", "_anc_aa"))
    old = {r.id: str(r.seq) for r in SeqIO.parse(aa_path, "fasta")}

    s = nt[a.name]
    o = old.get(a.name, "")
    print(f"{a.group.name}  {a.name}")
    print(f"  aligned nt len {len(s)}, ungapped {len(s.replace('-', ''))}, "
          f"first ATG at {s.replace('-', '').find('ATG')}")
    print(f"  old aa len {len(o)}")
    for m in ("ungapped", "aligned", "legacy"):
        try:
            n = nt_to_aa(s, mode=m)
        except Exception as e:  # noqa: BLE001
            print(f"  {m:<9} ERROR {e}")
            continue
        pref = n.startswith(o) or o.startswith(n)
        print(f"  {m:<9} len {len(n):<5} prefix-compatible={pref}")
    n = nt_to_aa(s, mode=a.mode)
    k = next((i for i in range(min(len(o), len(n))) if o[i] != n[i]), None)
    print(f"  first differing residue index: {k}")
    if k is not None:
        lo, hi = max(0, k - 12), k + 12
        print(f"    old[{lo}:{hi}] {o[lo:hi]}")
        print(f"    new[{lo}:{hi}] {n[lo:hi]}")
    print(f"  old head {o[:45]}")
    print(f"  new head {n[:45]}")


if __name__ == "__main__":
    main()
