#!/usr/bin/env python3
"""
Why do some new-root trees carry so few VOC-bundle leaves?

check_voc_roots_overlap.py found e.g. south_africa_beta with 120 B.1.351 leaves
by Pango label but only 3 matching Beta's defining bundle. Two candidate causes:
sequencing ambiguity (N runs -> X in the AA translation) at the bundle sites, or
a bundle definition that is stricter than the lineage.

Prints, per population, the distribution of bundle hits and of ambiguous
residues at the bundle positions, across leaves labelled as the target lineage.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

CASES = [
    ("south_africa_beta", "Beta", "001"),
    ("uk_alpha", "Alpha", "001"),
    ("uk_alpha", "Alpha", "003"),
    ("india_delta", "Delta", "001"),
    ("usa_epsilon", "Epsilon", "001"),
    ("usa_iota", "Iota", "001"),
    ("colombia_mu", "Mu", "001"),
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


def parse_mut(m: str) -> tuple[int, str]:
    i = 1
    while i < len(m) and m[i].isdigit():
        i += 1
    return int(m[1:i]), m[i:]


def main() -> None:
    panel = json.loads((ROOT / "benchmarks" / "voc_threat_panel.json").read_text())
    vocs = {v["id"]: v for v in panel["vocs"]}
    base = ROOT / "data" / "covid_voc_roots"

    for pop, voc_id, gid in CASES:
        fasta = base / pop / f"group_{gid}_anc_aa.fasta"
        meta = base / "raw" / f"{pop}_metadata.csv"
        if not fasta.exists() or not meta.exists():
            print(f"{pop} g{gid}: missing inputs, skipping")
            continue

        with meta.open() as fh:
            roles = {r["acc_base"]: r.get("role", "") for r in csv.DictReader(fh)}
        voc = vocs[voc_id]
        bundle = [parse_mut(m) for m in voc["defining_bundle"]]

        targets = [
            (n, s)
            for n, s in iter_fasta(fasta)
            if roles.get(n.split(",")[0]) == "target"
        ]
        print(f"\n=== {pop} g{gid} / {voc_id} ===")
        print(f"  bundle={voc['defining_bundle']} min_hits={voc.get('bundle_min_hits')}")
        print(f"  target-lineage leaves: {len(targets)}")
        if not targets:
            continue

        hits = Counter()
        ambig = Counter()
        per_site = Counter()
        per_site_ambig = Counter()
        for _, seq in targets:
            h = x = 0
            for pos, aa in bundle:
                c = seq[pos - 1] if pos <= len(seq) else "-"
                if c == aa:
                    h += 1
                    per_site[f"{pos}{aa}"] += 1
                if c in ("X", "-"):
                    x += 1
                    per_site_ambig[f"{pos}{aa}"] += 1
            hits[h] += 1
            ambig[x] += 1
        print(f"  bundle hits per leaf : {dict(sorted(hits.items()))}")
        print(f"  ambiguous sites/leaf : {dict(sorted(ambig.items()))}")
        print(f"  present per site     : {dict(per_site)}")
        print(f"  ambiguous per site   : {dict(per_site_ambig)}")


if __name__ == "__main__":
    main()
