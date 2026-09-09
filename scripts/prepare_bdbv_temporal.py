#!/usr/bin/env python3
"""
Temporal split for BDBV / pan-ebolavirus L window FASTA.

Default (tuned per DATA_AUDIT — sparse pre-2026 BDBV in NCBI L extract):
  train: BDBV year <= 2018
  val:   BDBV 2015 (optional; use --no-val if empty)
  test:  BDBV 2020-2026 (proxy; true 2026 holdout via Pathoplexus when available)

Output: data/bdbv_temporal/{train,val,test}/bdbvt{split}_group_NNN.fasta(+csv)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO

from scripts.bdbv_common import DATA_ROOT, format_fasta_header, parse_date, parse_fasta_header


def assign_split(
    species: str,
    year: int | None,
    train_end: int,
    val_year: int | None,
    test_start: int,
    test_end: int,
    pan: bool,
) -> str | None:
    if year is None:
        return None
    if species == "bdbv":
        if year <= train_end:
            return "train"
        if val_year is not None and year == val_year:
            return "val"
        if test_start <= year <= test_end:
            return "test"
        return None
    if pan:
        if test_start <= year <= test_end:
            return None  # non-BDBV never in test
        return "train"
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DATA_ROOT / "l_window" / "all.fasta"))
    ap.add_argument("--out-base", default="data/bdbv_temporal")
    ap.add_argument("--pan-ebolavirus", action="store_true")
    ap.add_argument("--train-end-year", type=int, default=2018)
    ap.add_argument("--val-year", type=int, default=2015)
    ap.add_argument("--no-val", action="store_true")
    ap.add_argument("--test-start-year", type=int, default=2020)
    ap.add_argument("--test-end-year", type=int, default=2026)
    ap.add_argument("--group-size", type=int, default=120)
    ap.add_argument("--min-group", type=int, default=5)
    args = ap.parse_args()

    inp = Path(args.input)
    if not inp.is_file():
        print(f"Missing {inp}")
        sys.exit(1)

    val_year = None if args.no_val else args.val_year
    pools: dict[str, list] = {"train": [], "val": [], "test": []}

    for rec in SeqIO.parse(inp, "fasta"):
        acc, date, country, species, source = parse_fasta_header(rec)
        _, _, year = parse_date(date)
        split = assign_split(
            species, year, args.train_end_year, val_year,
            args.test_start_year, args.test_end_year, args.pan_ebolavirus,
        )
        if split is None:
            continue
        pools[split].append((acc, date, country, species, source, str(rec.seq)))

    base = ROOT / args.out_base
    protocol = {
        "dataset": "bdbv_l",
        "split_type": "temporal",
        "pan_ebolavirus": args.pan_ebolavirus,
        "train_end_year": args.train_end_year,
        "val_year": val_year,
        "test_start_year": args.test_start_year,
        "test_end_year": args.test_end_year,
        "group_size": args.group_size,
        "counts": {k: len(v) for k, v in pools.items()},
    }
    base.mkdir(parents=True, exist_ok=True)
    (base / "SPLIT_PROTOCOL.json").write_text(json.dumps(protocol, indent=2) + "\n")

    prefixes = {"train": "bdbvttrain", "val": "bdbvtval", "test": "bdbvttest"}
    for split, recs in pools.items():
        if not recs:
            print(f"{split}: EMPTY")
            continue
        out_dir = base / split
        out_dir.mkdir(parents=True, exist_ok=True)
        for old in out_dir.glob(f"{prefixes[split]}_group_*"):
            old.unlink()
        pool = out_dir / f"{prefixes[split]}.fasta"
        recs.sort(key=lambda x: x[1])
        with open(pool, "w") as f:
            for acc, date, country, species, source, seq in recs:
                f.write(format_fasta_header(acc, date, country, species, source) + f"\n{seq}\n")
        print(f"{split}: {len(recs)} seqs -> grouping")
        for gi in range(0, len(recs), args.group_size):
            chunk = recs[gi : gi + args.group_size]
            if len(chunk) < args.min_group and split != "test":
                continue
            gnum = (gi // args.group_size) + 1
            gf = out_dir / f"{prefixes[split]}_group_{gnum:03d}.fasta"
            gcsv = out_dir / f"{prefixes[split]}_group_{gnum:03d}.csv"
            with open(gf, "w") as ff, open(gcsv, "w", newline="") as cf:
                import csv

                w = csv.writer(cf)
                w.writerow(["name", "date", "country", "species"])
                for acc, date, country, species, source, seq in chunk:
                    ff.write(format_fasta_header(acc, date, country, species, source) + f"\n{seq}\n")
                    w.writerow([acc, date, country, species])
            print(f"  group {gnum}: {len(chunk)} seqs")

    print(f"Done -> {base}")


if __name__ == "__main__":
    main()
