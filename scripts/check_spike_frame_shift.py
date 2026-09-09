#!/usr/bin/env python3
"""
Test whether Spike deletions shift downstream reference coordinates.

covid_extract_spike.py slices the aligned genome at Spike's reference
coordinates and then strips gaps (`.replace("-", "")`). Stripping gaps is what
breaks the frame: a lineage carrying a Spike deletion ends up shorter, so every
residue downstream of the deletion sits at a lower index than its reference
position. Mutation calls that index by reference position then silently miss.

Alpha (dH69/V70 + dY144, 3 residues) and Beta (d242-244, 3 residues) both delete
upstream of their defining sites, which would explain why only a handful of
B.1.1.7/B.1.351-labelled leaves appear to carry N501Y or K417N.

Prints, per case, how many target leaves carry each defining mutation at the
reference position versus at positions shifted by 1..6.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

CASES = [
    ("south_africa_beta", "Beta", "001"),
    ("uk_alpha", "Alpha", "001"),
    ("usa_epsilon", "Epsilon", "001"),
]


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


def parse_mut(m: str) -> tuple[str, int, str]:
    i = 1
    while i < len(m) and m[i].isdigit():
        i += 1
    return m[0], int(m[1:i]), m[i:]


def main() -> None:
    panel = json.loads((ROOT / "benchmarks" / "voc_threat_panel.json").read_text())
    vocs = {v["id"]: v for v in panel["vocs"]}
    base = ROOT / "data" / "covid_voc_roots"
    wt_path = ROOT / "data" / "covid" / "wt.txt"
    wt = wt_path.read_text().strip() if wt_path.exists() else ""

    for pop, voc_id, gid in CASES:
        fasta = base / pop / f"group_{gid}_anc_aa.fasta"
        meta = base / "raw" / f"{pop}_metadata.csv"
        if not fasta.exists():
            print(f"{pop}: missing {fasta}")
            continue
        with meta.open() as fh:
            roles = {r["acc_base"]: r.get("role", "") for r in csv.DictReader(fh)}
        targets = [
            s for n, s in iter_fasta(fasta) if roles.get(n.split(",")[0]) == "target"
        ]
        voc = vocs[voc_id]
        print(f"\n=== {pop} / {voc_id} ({len(targets)} target leaves) ===")
        lens = sorted({len(s) for s in targets})
        print(f"  distinct leaf lengths: {lens[:8]}{' …' if len(lens) > 8 else ''}")
        if wt:
            print(f"  reference (wt.txt) length: {len(wt)}")

        for mut in voc["defining_bundle"]:
            ref_aa, pos, alt = parse_mut(mut)
            row = []
            for shift in range(0, 7):
                idx = pos - 1 - shift
                n = sum(1 for s in targets if 0 <= idx < len(s) and s[idx] == alt)
                row.append(f"-{shift}:{n}" if shift else f"ref:{n}")
            print(f"  {mut:<8} " + "  ".join(row))


if __name__ == "__main__":
    main()
