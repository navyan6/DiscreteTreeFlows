#!/usr/bin/env python3
"""
Outbreak-aware pan-filovirus L split for TreeSBM.

One tree = one outbreak (single species). Large outbreaks split into date-ordered
sub-chunks of --group-size within the same outbreak_id.

Track A (default):
  train  = all outbreaks except val + test
  val    = one held-out EBOV outbreak (--val-outbreak auto)
  test   = BDBV 2026 only (outbreak bdbv_2026 or BDBV year >= 2026)

Track B (--track-b):
  test   = --test-outbreak (e.g. ebov_wa_2013_2016) for rich sanity eval

Output:
  data/filo_l/{train,val,test}/filo_{split}_group_NNN.fasta(+csv)
  data/filo_l/SPLIT_PROTOCOL.json

Prereq:
  python scripts/bdbv_extract_l.py
  python scripts/bdbv_define_l_window.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO

from scripts.bdbv_common import (
    DATA_ROOT,
    GOLDEN_TEST_OUTBREAK,
    GOLDEN_TEST_OUTBREAKS,
    format_fasta_header,
    infer_outbreak_id,
    load_manifest,
    parse_date,
    parse_fasta_header,
    slug_outbreak_id,
)

DEFAULT_OUT_BASE = "data/filo_l"
DEFAULT_VAL_OUTBREAK = "auto"
DEFAULT_TEST_OUTBREAK_TRACK_A = GOLDEN_TEST_OUTBREAK
PREFIXES = {"train": "filo_train", "val": "filo_val", "test": "filo_test"}


def _load_records(input_fasta: Path) -> list[dict]:
    manifest = {m["accession"].split(".")[0]: m for m in load_manifest()}
    records: list[dict] = []
    for rec in SeqIO.parse(input_fasta, "fasta"):
        acc, date, country, species, source = parse_fasta_header(rec)
        acc = acc.split(".")[0]
        _, _, year = parse_date(date)
        meta = manifest.get(acc, {})
        if meta.get("golden_test") or meta.get("outbreak_id") in GOLDEN_TEST_OUTBREAKS:
            outbreak_id = GOLDEN_TEST_OUTBREAK
        else:
            outbreak_id = meta.get("outbreak_id") or infer_outbreak_id(
                {
                    "accession": acc,
                    "species": species,
                    "country": country,
                    "year": year,
                    "outbreak": meta.get("outbreak", ""),
                }
            )
        records.append(
            {
                "acc": acc,
                "date": date,
                "country": country,
                "species": species,
                "source": source,
                "year": year,
                "outbreak_id": outbreak_id,
                "seq": str(rec.seq),
            }
        )
    return records


def _pick_val_outbreak(
    by_outbreak: dict[str, list[dict]],
    reserved: set[str],
    target_n: int = 60,
) -> str | None:
    candidates: list[tuple[int, str]] = []
    for oid, recs in by_outbreak.items():
        if oid in reserved:
            continue
        if recs[0]["species"] != "ebov":
            continue
        n = len(recs)
        if n < 15:
            continue
        candidates.append((n, oid))
    if not candidates:
        return None
    candidates.sort(key=lambda x: (abs(x[0] - target_n), -x[0]))
    return candidates[0][1]


def _assign_outbreak_splits(
    outbreak_ids: set[str],
    track_b: bool,
    test_outbreak: str,
    val_outbreak: str,
    by_outbreak: dict[str, list[dict]],
) -> dict[str, str]:
    test_oid = slug_outbreak_id(test_outbreak)
    val_oid = slug_outbreak_id(val_outbreak) if val_outbreak != "auto" else None

    if track_b:
        test_set = {test_oid}
    else:
        test_set = {
            oid
            for oid in outbreak_ids
            if oid == DEFAULT_TEST_OUTBREAK_TRACK_A
            or (
                by_outbreak[oid][0]["species"] == "bdbv"
                and by_outbreak[oid][0].get("year") is not None
                and by_outbreak[oid][0]["year"] >= 2026
            )
        }

    # Gold Pathoplexus / BDBV 2026 band — eval only, never train or val.
    test_set |= GOLDEN_TEST_OUTBREAKS & outbreak_ids

    if val_oid is None:
        val_oid = _pick_val_outbreak(by_outbreak, reserved=test_set)
    val_set = {val_oid} if val_oid else set()
    val_set -= test_set

    split_map: dict[str, str] = {}
    for oid in outbreak_ids:
        if oid in test_set:
            split_map[oid] = "test"
        elif oid in val_set:
            split_map[oid] = "val"
        else:
            split_map[oid] = "train"
    return split_map


def _chunk_outbreak(recs: list[dict], group_size: int, min_group: int) -> list[list[dict]]:
    recs = sorted(recs, key=lambda r: r["date"])
    if len(recs) <= group_size:
        return [recs] if len(recs) >= min_group else []

    chunks: list[list[dict]] = []
    for i in range(0, len(recs), group_size):
        chunk = recs[i : i + group_size]
        if len(chunk) >= min_group or (chunks and i + group_size >= len(recs)):
            chunks.append(chunk)
        elif not chunks and len(chunk) >= max(3, min_group // 2):
            chunks.append(chunk)
    return chunks


def _write_groups(
    pools: dict[str, list[list[dict]]],
    out_base: Path,
    min_group_test: int,
) -> dict[str, int]:
    stats: dict[str, int] = {"train": 0, "val": 0, "test": 0}
    for split, group_lists in pools.items():
        out_dir = out_base / split
        out_dir.mkdir(parents=True, exist_ok=True)
        for old in out_dir.glob(f"{PREFIXES[split]}_group_*"):
            old.unlink()
        gnum = 0
        prefix = PREFIXES[split]
        min_g = 3 if split == "test" else min_group_test
        for chunk in group_lists:
            if len(chunk) < min_g and split != "test":
                continue
            if len(chunk) < 2 and split == "test":
                continue
            gnum += 1
            gf = out_dir / f"{prefix}_group_{gnum:03d}.fasta"
            gcsv = out_dir / f"{prefix}_group_{gnum:03d}.csv"
            with open(gf, "w") as ff, open(gcsv, "w", newline="") as cf:
                w = csv.writer(cf)
                w.writerow(["name", "date", "country", "species", "outbreak_id"])
                for r in chunk:
                    ff.write(
                        format_fasta_header(
                            r["acc"], r["date"], r["country"], r["species"], r["source"]
                        )
                        + f"\n{r['seq']}\n"
                    )
                    w.writerow(
                        [r["acc"], r["date"], r["country"], r["species"], r["outbreak_id"]]
                    )
            stats[split] += len(chunk)
            print(
                f"  {split} group {gnum:03d}: {len(chunk)} seqs "
                f"outbreak={chunk[0]['outbreak_id']} species={chunk[0]['species']}"
            )
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DATA_ROOT / "l_window" / "all.fasta"))
    ap.add_argument("--out-base", default=DEFAULT_OUT_BASE)
    ap.add_argument("--group-size", type=int, default=80)
    ap.add_argument("--min-group", type=int, default=5)
    ap.add_argument("--val-outbreak", default=DEFAULT_VAL_OUTBREAK)
    ap.add_argument(
        "--test-outbreak",
        default=DEFAULT_TEST_OUTBREAK_TRACK_A,
        help="Track B only; Track A always uses bdbv_2026",
    )
    ap.add_argument(
        "--track-b",
        action="store_true",
        help="Rich holdout: test=--test-outbreak instead of BDBV 2026",
    )
    args = ap.parse_args()

    inp = Path(args.input)
    if not inp.is_file():
        print(f"Missing {inp}")
        sys.exit(1)

    records = _load_records(inp)
    by_outbreak: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        by_outbreak[r["outbreak_id"]].append(r)

    # Single-species guard: split mixed outbreaks by species suffix
    cleaned: dict[str, list[dict]] = {}
    for oid, recs in by_outbreak.items():
        species_set = {r["species"] for r in recs}
        if len(species_set) == 1:
            cleaned[oid] = recs
            continue
        for sp in sorted(species_set):
            sub = [r for r in recs if r["species"] == sp]
            cleaned[slug_outbreak_id(f"{oid}_{sp}")] = sub

    outbreak_ids = set(cleaned.keys())
    split_map = _assign_outbreak_splits(
        outbreak_ids,
        args.track_b,
        args.test_outbreak,
        args.val_outbreak,
        cleaned,
    )

    pools: dict[str, list[list[dict]]] = {"train": [], "val": [], "test": []}
    outbreak_hist: dict[str, dict] = {}
    for oid, recs in sorted(cleaned.items()):
        split = split_map.get(oid, "train")
        if oid in GOLDEN_TEST_OUTBREAKS:
            split = "test"
        elif split == "train" and any(r["outbreak_id"] in GOLDEN_TEST_OUTBREAKS for r in recs):
            split = "test"
        chunks = _chunk_outbreak(recs, args.group_size, args.min_group)
        pools[split].extend(chunks)
        outbreak_hist[oid] = {
            "split": split,
            "n_seqs": len(recs),
            "n_groups": len(chunks),
            "species": recs[0]["species"],
        }

    base = ROOT / args.out_base
    base.mkdir(parents=True, exist_ok=True)

    print(f"Writing groups -> {base}")
    counts = _write_groups(pools, base, args.min_group)

    val_outbreaks = [o for o, s in split_map.items() if s == "val"]
    test_outbreaks = [o for o, s in split_map.items() if s == "test"]
    protocol = {
        "dataset": "filo_l",
        "split_type": "outbreak",
        "track": "B" if args.track_b else "A",
        "group_size": args.group_size,
        "min_group": args.min_group,
        "val_outbreak": val_outbreaks[0] if val_outbreaks else None,
        "test_outbreaks": test_outbreaks,
        "tree_semantics": "one FastTree per outbreak_id (single species); "
        "large outbreaks split by date within outbreak",
        "counts": counts,
        "n_outbreaks": len(outbreak_hist),
        "outbreaks": outbreak_hist,
    }
    (base / "SPLIT_PROTOCOL.json").write_text(json.dumps(protocol, indent=2) + "\n")

    for split in ("train", "val", "test"):
        ng = len(list((base / split).glob("*_group_*.fasta"))) if (base / split).is_dir() else 0
        print(f"{split}: {counts[split]} seqs in {ng} groups")
    if counts["val"] == 0:
        print("WARNING: val band empty — training will fall back to test for early stopping")
    if counts["test"] == 0:
        print("WARNING: test band empty — add Pathoplexus 2026 or lower test year threshold")


if __name__ == "__main__":
    main()
