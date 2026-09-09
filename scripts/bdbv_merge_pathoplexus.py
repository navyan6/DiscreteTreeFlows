#!/usr/bin/env python3
"""Merge Pathoplexus FASTA into manifest (source=pathoplexus)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO

from scripts.bdbv_common import (
    DATA_ROOT,
    dedupe_manifest,
    load_manifest,
    normalize_country,
    parse_date,
    save_manifest,
    write_data_audit,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", default=str(DATA_ROOT / "raw" / "pathoplexus_2026"))
    args = ap.parse_args()
    indir = Path(args.input_dir)
    new = []
    for fp in sorted(indir.glob("**/*")):
        if fp.suffix.lower() not in (".fasta", ".fa", ".gb", ".gbk"):
            continue
        fmt = "genbank" if fp.suffix.lower() in (".gb", ".gbk") else "fasta"
        for rec in SeqIO.parse(fp, fmt):
            acc = rec.id.split(".")[0]
            desc = rec.description or ""
            parts = desc.split(",")
            date_raw = parts[1] if len(parts) > 1 else "2026-01-01"
            country = normalize_country(parts[2] if len(parts) > 2 else "uganda")
            augur_date, _, year = parse_date(date_raw)
            new.append(
                {
                    "accession": acc,
                    "species": "bdbv",
                    "date": augur_date,
                    "year": year,
                    "country": country,
                    "source": "pathoplexus",
                    "length": len(rec.seq),
                    "host": "Homo sapiens",
                    "fasta": str(fp.relative_to(ROOT)),
                }
            )
    merged = dedupe_manifest(load_manifest() + new)
    save_manifest(merged)
    write_data_audit(merged)
    print(f"Added {len(new)} pathoplexus records; manifest size {len(merged)}")


if __name__ == "__main__":
    main()
