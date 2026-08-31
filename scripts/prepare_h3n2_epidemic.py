#!/usr/bin/env python3
"""
Epidemic-aware H3N2 HA split: one tree per NH flu season (+ date chunks if huge).

Temporal cutoffs match prepare_h3n2_temporal.py defaults (train≤2022 / val2023 / test2024).
Groups by Northern-hemisphere season (Oct–Sep), not arbitrary 400-seq date batches.

Output: data/h3n2_epidemic/{split}/h3n2etrain_group_NNN.fasta(+csv)
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
    write_fasta_csv_groups,
    write_split_protocol,
)
from scripts.prepare_h3n2_temporal import SOURCE_WINDOWS, _id, _year

PREFIXES = {"train": "h3n2etrain", "val": "h3n2eval", "test": "h3n2etest"}
POOL_NAMES = {"train": "h3n2train", "val": "h3n2val", "test": "h3n2test"}


def _load_from_pools(pool_base: Path) -> tuple[dict[str, dict[str, list]], dict]:
    """Regroup existing temporal split pools by NH season (lab pipeline outputs)."""
    pools: dict[str, dict[str, list]] = {
        "train": defaultdict(list),
        "val": defaultdict(list),
        "test": defaultdict(list),
    }
    season_hist: dict[str, dict[str, int]] = {s: {} for s in pools}
    for split, pool_name in POOL_NAMES.items():
        pool = pool_base / split / f"{pool_name}.fasta"
        if not pool.is_file():
            print(f"WARNING: missing {pool}")
            continue
        for rec in SeqIO.parse(pool, "fasta"):
            rid, date_str, _, _, season = parse_flu_header(rec.description)
            if not season:
                try:
                    season = f"unknown_{int(date_str[:4])}"
                except ValueError:
                    season = "unknown"
            sort_key = date_str
            row = (rid, date_str, sort_key, str(rec.seq), season)
            pools[split][season].append(row)
            season_hist[split][season] = season_hist[split].get(season, 0) + 1
        print(f"  {split}: {sum(season_hist[split].values())} seqs from {pool.name}")
    return pools, {"season_hist": season_hist, "source": "pool_base"}


def _load_pools(
    train_end: int,
    val_year: int | None,
    test_start: int,
    test_end: int,
    min_year: int | None,
) -> tuple[dict[str, dict[str, list]], dict]:
    pools: dict[str, dict[str, list]] = {
        "train": defaultdict(list),
        "val": defaultdict(list),
        "test": defaultdict(list),
    }
    seen: set[str] = set()
    season_hist: dict[str, dict[str, int]] = {s: {} for s in pools}
    window_map: list[dict] = []

    for wp, y_lo, y_hi in SOURCE_WINDOWS:
        path = ROOT / wp
        entry = {"path": wp, "exists": path.is_file()}
        if not path.is_file():
            window_map.append(entry)
            continue
        for rec in SeqIO.parse(path, "fasta"):
            y = _year(rec)
            if y is None:
                continue
            if min_year is not None and y < min_year:
                continue
            split = assign_temporal_split(y, train_end, val_year, test_start, test_end)
            if split is None:
                continue
            rid = _id(rec)
            if rid in seen:
                continue
            seen.add(rid)
            _, date_str, _, _, season = parse_flu_header(rec.description)
            if not season:
                season = f"unknown_{y}"
            sort_key = date_str
            row = (rid, date_str, sort_key, str(rec.seq), season)
            pools[split][season].append(row)
            season_hist[split][season] = season_hist[split].get(season, 0) + 1
        window_map.append(entry)

    return pools, {"windows": window_map, "season_hist": season_hist}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-base", default="data/h3n2_epidemic")
    ap.add_argument("--group-size", type=int, default=400)
    ap.add_argument("--min-group", type=int, default=40)
    ap.add_argument("--min-year", type=int, default=2014)
    ap.add_argument("--train-end-year", type=int, default=2022)
    ap.add_argument("--val-year", type=int, default=2023)
    ap.add_argument("--no-val", action="store_true")
    ap.add_argument("--test-start-year", type=int, default=2024)
    ap.add_argument("--test-end-year", type=int, default=2024)
    ap.add_argument(
        "--pool-base",
        default="",
        help="Use existing temporal pools e.g. lab .../h3n2 (train/h3n2train.fasta)",
    )
    args = ap.parse_args()

    val_year = None if args.no_val else args.val_year
    base = ROOT / args.out_base

    if args.pool_base:
        pools, meta = _load_from_pools(ROOT / args.pool_base)
    else:
        pools, meta = _load_pools(
            args.train_end_year,
            val_year,
            args.test_start_year,
            args.test_end_year,
            args.min_year,
        )

    counts = {}
    n_groups = {}
    for split, prefix in PREFIXES.items():
        out_dir = base / split
        for old in out_dir.glob(f"{prefix}_group_*"):
            old.unlink()
        g, total = 1, 0
        for season in sorted(pools[split].keys()):
            chunks = chunk_records(
                pools[split][season], args.group_size, args.min_group, sort_key_idx=2
            )
            if not chunks:
                continue
            g_end, n = write_fasta_csv_groups(
                chunks,
                out_dir,
                prefix,
                g,
                extra_csv_cols=["season"],
                extra_row_fn=lambda r: [r[0], r[1], r[4]],
            )
            print(f"  [{split}] season {season}: {len(pools[split][season])} seqs -> "
                  f"groups {g:03d}-{g_end - 1:03d}")
            g = g_end
            total += n
        counts[split] = total
        n_groups[split] = g - 1
        print(f"=== {split}: {total} seqs, {g - 1} groups ===\n")

    protocol = {
        "dataset": "h3n2",
        "split_type": "epidemic",
        "out_base": args.out_base,
        "train_end_year": args.train_end_year,
        "val_year": val_year,
        "test_start_year": args.test_start_year,
        "test_end_year": args.test_end_year,
        "min_year": args.min_year,
        "group_size": args.group_size,
        "tree_semantics": "one tree per NH flu season (Oct(Y-1)–Sep(Y)); "
        "large seasons split by date within season",
        "season_hist": meta["season_hist"],
        "counts": counts,
        "n_groups": n_groups,
        "checkpoint_target": "checkpoints/h3n2_v4_epidemic_mutrec",
    }
    write_split_protocol(base, protocol)
    print(f"Wrote {base / 'SPLIT_PROTOCOL.json'}")


if __name__ == "__main__":
    main()
