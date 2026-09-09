#!/usr/bin/env python3
"""
Geographic split for BDBV L window: East Africa holdout vs rest.

Test = East Africa (Uganda, Sudan, DRC); val = 10% of train units; train = rest.
Pan mode: only BDBV in val/test bands; other species -> train only.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO

from scripts.bdbv_common import DATA_ROOT, is_east_africa, parse_date, parse_fasta_header, format_fasta_header


def greedy_split(counts: dict[str, int]) -> dict[str, str]:
    total = sum(counts.values())
    target = {"train": 0.80 * total, "val": 0.10 * total, "test": 0.10 * total}
    fill = {"train": 0.0, "val": 0.0, "test": 0.0}
    assign = {}
    for unit, c in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        split = max(("train", "val", "test"), key=lambda s: target[s] - fill[s])
        assign[unit] = split
        fill[split] += c
    return assign


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DATA_ROOT / "l_window" / "all.fasta"))
    ap.add_argument("--out-base", default="data/bdbv_geo")
    ap.add_argument("--pan-ebolavirus", action="store_true")
    ap.add_argument("--group-size", type=int, default=120)
    ap.add_argument("--min-group", type=int, default=20)
    args = ap.parse_args()

    inp = Path(args.input)
    if not inp.is_file():
        sys.exit(f"Missing {inp}")

    by_unit: dict[str, list] = {}
    for rec in SeqIO.parse(inp, "fasta"):
        acc, date, country, species, source = parse_fasta_header(rec)
        if species != "bdbv" and not args.pan_ebolavirus:
            continue
        if species != "bdbv" and args.pan_ebolavirus:
            unit = f"{species}_global"
            split_band = "train"
        else:
            unit = country or "unknown"
            split_band = "test" if is_east_africa(country) else "train_pool"
        by_unit.setdefault(unit, []).append(
            (acc, date, country, species, source, str(rec.seq), split_band)
        )

    # assign val from train_pool units
    train_pool_units = {
        u: len([x for x in recs if x[6] == "train_pool"])
        for u, recs in by_unit.items()
        if any(x[6] == "train_pool" for x in recs)
    }
    unit_split = greedy_split(train_pool_units)
    pools: dict[str, list] = {"train": [], "val": [], "test": []}

    for unit, recs in by_unit.items():
        for acc, date, country, species, source, seq, band in recs:
            row = (acc, date, country, species, source, seq)
            if band == "test":
                pools["test"].append(row)
            elif band == "train_pool":
                sp = unit_split.get(unit, "train")
                pools[sp].append(row)
            else:
                pools["train"].append(row)

    base = ROOT / args.out_base
    base.mkdir(parents=True, exist_ok=True)
    (base / "SPLIT_PROTOCOL.json").write_text(
        json.dumps(
            {
                "split_type": "geographic",
                "pan_ebolavirus": args.pan_ebolavirus,
                "test_region": "east_africa",
                "counts": {k: len(v) for k, v in pools.items()},
            },
            indent=2,
        )
        + "\n"
    )

    prefixes = {"train": "bdbvgtrain", "val": "bdbvgval", "test": "bdbvgtest"}
    for split, recs in pools.items():
        if not recs:
            continue
        out_dir = base / split
        out_dir.mkdir(parents=True, exist_ok=True)
        recs.sort(key=lambda x: x[1])
        g = 0
        for i in range(0, len(recs), args.group_size):
            chunk = recs[i : i + args.group_size]
            if len(chunk) < args.min_group:
                continue
            g += 1
            fasta = out_dir / f"{prefixes[split]}_group_{g:03d}.fasta"
            csvp = out_dir / f"{prefixes[split]}_group_{g:03d}.csv"
            with open(fasta, "w") as ff, open(csvp, "w", newline="") as cf:
                w = csv.writer(cf)
                w.writerow(["name", "date"])
                for acc, date, country, species, source, seq in chunk:
                    ff.write(format_fasta_header(acc, date, country, species, source) + f"\n{seq}\n")
                    w.writerow([acc, date])
        print(f"{split}: {len(recs)} seqs -> {g} groups")


if __name__ == "__main__":
    main()
