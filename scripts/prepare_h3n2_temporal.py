#!/usr/bin/env python3
"""
Prepare the temporally-split H3N2 HA dataset for the TreeSBM pipeline.

Splits (confirmed design):
  train = 2014-2022 windows   (2023 dropped to avoid val leakage)
  val   = 2023
  test  = 2024 (early + full)

For each split we concatenate the source windows, deduplicate by EPI_ISL id, and
enforce fully disjoint splits (val \\ train, test \\ (train ∪ val)) so there is no
sequence leakage. Each split pool is then date-ordered and cut into groups of
`--group-size` (default 400) via split_fasta_by_date, written to its own dir so
the flat group_NNN numbering never collides with the existing multi-subtype data
in data/train/.

Output layout (gitignore this):
  data/h3n2/train/h3n2train_group_NNN.fasta (+ .csv)
  data/h3n2/val/h3n2val_group_NNN.fasta     (+ .csv)
  data/h3n2/test/h3n2test_group_NNN.fasta   (+ .csv)

Downstream:
  run_all_groups.py --data-dir data/h3n2/train --prefix h3n2train   (+ val, test)
  precompute_plm.py / precompute_ref_rates.py --data data/h3n2/{split}
  train.py --data data/h3n2/train --val-data data/h3n2/val --test-data data/h3n2/test
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO
from scripts.split_fasta_by_date import split_fasta_by_date

# Source windows (relative to repo root)
TRAIN_WINDOWS = [
    "data/train/h3n2_train/2014_2016_h3n2_ha.fasta",
    "data/train/h3n2_train/2016_2017_h3n2_ha.fasta",
    "data/train/h3n2_train/2017_2018_h3n2_ha.fasta",
    "data/train/h3n2_train/2018_2019_h3n2_ha.fasta",
    "data/train/h3n2_train/2019_2020_h3n2_ha.fasta",
    "data/train/h3n2_train/2020_2022_h3n2_ha.fasta",
    "data/train/h3n2_train/2022_h3n2_ha.fasta",
]
VAL_WINDOWS = ["data/validate/h3n2_val/2023_h3n2_ha.fasta"]
TEST_WINDOWS = [
    "data/test/h3n2_ha_test/h3n2_ha_early_2024.fasta",
    "data/test/h3n2_ha_test/h3n2_ha_2024.fasta",
]


def _id(rec) -> str:
    # headers: >EPI_ISL_XXXXXX,YYYY-MM-DD
    return rec.description.split(",", 1)[0].strip()


def load_pool(window_paths, exclude_ids: set[str]) -> tuple[list, set[str]]:
    """Concatenate windows, dedup by id, drop ids in `exclude_ids`. Returns (records, ids)."""
    seen: set[str] = set()
    records = []
    dropped_dup = dropped_excl = 0
    for wp in window_paths:
        path = ROOT / wp
        if not path.exists():
            print(f"  WARNING: missing {wp}")
            continue
        for rec in SeqIO.parse(path, "fasta"):
            rid = _id(rec)
            if rid in exclude_ids:
                dropped_excl += 1
                continue
            if rid in seen:
                dropped_dup += 1
                continue
            seen.add(rid)
            records.append(rec)
    print(f"  kept {len(records)} | dropped {dropped_dup} intra-split dups, "
          f"{dropped_excl} cross-split leaks")
    return records, seen


def write_pool(records, out_fasta: Path):
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fasta, "w") as f:
        for rec in records:
            # preserve the ">id,date" header split_fasta_by_date expects
            f.write(f">{rec.description}\n{str(rec.seq)}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-base", default="data/h3n2")
    ap.add_argument("--group-size", type=int, default=400)
    args = ap.parse_args()

    base = ROOT / args.out_base

    print("[train] 2014-2022 windows")
    train_recs, train_ids = load_pool(TRAIN_WINDOWS, exclude_ids=set())
    print("[val] 2023 (minus train ids)")
    val_recs, val_ids = load_pool(VAL_WINDOWS, exclude_ids=train_ids)
    print("[test] 2024 (minus train+val ids)")
    test_recs, _ = load_pool(TEST_WINDOWS, exclude_ids=train_ids | val_ids)

    splits = {
        "train": ("h3n2train", train_recs),
        "val":   ("h3n2val",   val_recs),
        "test":  ("h3n2test",  test_recs),
    }
    for split, (prefix, recs) in splits.items():
        out_dir = base / split
        out_dir.mkdir(parents=True, exist_ok=True)
        pool = out_dir / f"{prefix}.fasta"
        write_pool(recs, pool)
        print(f"\n=== splitting {split}: {len(recs)} seqs -> groups of {args.group_size} ===")
        split_fasta_by_date(str(pool), args.group_size, str(out_dir))

    print("\nDone. Next: run_all_groups.py --data-dir on each split dir.")


if __name__ == "__main__":
    main()
