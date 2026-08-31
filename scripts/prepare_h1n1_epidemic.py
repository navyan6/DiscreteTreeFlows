#!/usr/bin/env python3
"""
Epidemic-aware H1N1 HA split: location × NH flu season trees.

Temporal cutoffs match prepare_h1n1_temporal.py (train≤2023 / val2024 / test2025).
Each tree = (location_or_country, flu season) — coherent local season epidemic.

Output: data/h1n1_epidemic/{split}/h1n1etrain_group_NNN.fasta(+csv)
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO

from scripts.epidemic_split_common import (
    assign_temporal_split,
    chunk_records,
    parse_flu_header,
    slug_label,
    write_fasta_csv_groups,
    write_split_protocol,
)
from scripts.prepare_h1n1_geo import RAW_DIR, normalize_location, parse_date, parse_header

LOC_MIN = 100
GROUP_SIZE = 300
MIN_GROUP = 50
PREFIXES = {"train": "h1n1etrain", "val": "h1n1eval", "test": "h1n1etest"}
POOL_NAMES = {"train": "h1n1ttrain", "val": "h1n1tval", "test": "h1n1ttest"}


def _record_header(rec) -> str:
    desc = (rec.description or "").strip()
    if desc.startswith("|"):
        return f"{rec.id}{desc}"
    if "|" in desc:
        return desc
    return f"{rec.id} {desc}".strip()


def _load_from_pools(pool_base: Path) -> dict[str, dict[str, list]]:
    pools: dict[str, dict[str, list]] = {
        "train": defaultdict(list),
        "val": defaultdict(list),
        "test": defaultdict(list),
    }
    for split, pool_name in POOL_NAMES.items():
        pool = pool_base / split / f"{pool_name}.fasta"
        if not pool.is_file():
            print(f"WARNING: missing {pool}")
            continue
        for rec in SeqIO.parse(pool, "fasta"):
            acc = rec.id.split()[0].split(",")[0]
            date_str = rec.description.strip() if rec.description else "2000-01-01"
            if "," in rec.id:
                date_str = rec.id.split(",", 1)[1].strip()
            _, _, _, _, season = parse_flu_header(f"{acc},{date_str}")
            if not season:
                season = f"unknown_{date_str[:4]}"
            pools[split][season].append(
                (acc, date_str, date_str, str(rec.seq), "season", season, "")
            )
        print(f"  {split}: {sum(len(v) for v in pools[split].values())} seqs from {pool.name}")
    return pools


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-base", default="data/h1n1_epidemic")
    ap.add_argument("--loc-min", type=int, default=LOC_MIN)
    ap.add_argument("--group-size", type=int, default=GROUP_SIZE)
    ap.add_argument("--min-group", type=int, default=MIN_GROUP)
    ap.add_argument("--min-year", type=int, default=2009)
    ap.add_argument("--train-end-year", type=int, default=2023)
    ap.add_argument("--val-year", type=int, default=2024)
    ap.add_argument("--no-val", action="store_true")
    ap.add_argument("--test-start-year", type=int, default=2025)
    ap.add_argument("--test-end-year", type=int, default=2025)
    ap.add_argument(
        "--pool-base",
        default="",
        help="Use lab temporal pools e.g. .../h1n1_temporal",
    )
    args = ap.parse_args()

    val_year = None if args.no_val else args.val_year
    base = ROOT / args.out_base

    if args.pool_base:
        pools = _load_from_pools(ROOT / args.pool_base)
    else:
        raw_dir = ROOT / RAW_DIR
        raw_files = sorted(raw_dir.glob("*_h1n1.fasta"))
        if not raw_files:
            print(f"No FASTA in {raw_dir}")
            sys.exit(1)

        loc_counts: dict[str, dict[str, int]] = {
            "train": defaultdict(int),
            "val": defaultdict(int),
            "test": defaultdict(int),
        }
        raw_rows: list[tuple] = []
        seen: set[str] = set()

        for fp in raw_files:
            for rec in SeqIO.parse(fp, "fasta"):
                acc, location, date_raw, country = parse_header(_record_header(rec))
                if acc in seen:
                    continue
                try:
                    date_str, sort_key = parse_date(date_raw)
                except (ValueError, IndexError):
                    continue
                y = sort_key[0]
                if y < args.min_year:
                    continue
                split = assign_temporal_split(
                    y, args.train_end_year, val_year, args.test_start_year, args.test_end_year
                )
                if split is None:
                    continue
                seen.add(acc)
                loc = normalize_location(location)
                unit_key = loc if loc else country
                loc_counts[split][unit_key] += 1
                _, _, _, _, season = parse_flu_header(f"{acc},{date_raw}")
                if not season:
                    season = f"unknown_{y}"
                raw_rows.append((split, acc, date_str, sort_key, str(rec.seq), unit_key, country, season))

        pools = {
            "train": defaultdict(list),
            "val": defaultdict(list),
            "test": defaultdict(list),
        }
        for split, acc, date_str, sort_key, seq, unit_key, country, season in raw_rows:
            unit = unit_key if loc_counts[split][unit_key] >= args.loc_min else country
            ekey = slug_label(f"{unit}_{season}")
            pools[split][ekey].append((acc, date_str, sort_key, seq, unit, season, country))

    counts = {}
    n_groups = {}
    for split, prefix in PREFIXES.items():
        out_dir = base / split
        for old in out_dir.glob(f"{prefix}_group_*"):
            old.unlink()
        g, total = 1, 0
        for ekey in sorted(pools[split].keys()):
            chunks = chunk_records(
                pools[split][ekey], args.group_size, args.min_group, sort_key_idx=2
            )
            if not chunks:
                continue
            g_end, n = write_fasta_csv_groups(
                chunks,
                out_dir,
                prefix,
                g,
                extra_csv_cols=["unit", "season", "country"],
                extra_row_fn=lambda r: [r[0], r[1], r[4], r[5], r[6]],
            )
            g = g_end
            total += n
        counts[split] = total
        n_groups[split] = g - 1
        print(f"=== {split}: {total} seqs, {g - 1} groups ===\n")

    protocol = {
        "dataset": "h1n1",
        "split_type": "epidemic",
        "out_base": args.out_base,
        "train_end_year": args.train_end_year,
        "val_year": val_year,
        "test_start_year": args.test_start_year,
        "test_end_year": args.test_end_year,
        "tree_semantics": "one tree per (location|country, NH flu season)",
        "counts": counts,
        "n_groups": n_groups,
        "checkpoint_target": "checkpoints/h1n1_v3_epidemic_mutrec",
    }
    write_split_protocol(base, protocol)
    print(f"Wrote {base / 'SPLIT_PROTOCOL.json'}")


if __name__ == "__main__":
    main()
