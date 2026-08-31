#!/usr/bin/env python3
"""Audit filovirus outbreak split; gate training submit."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO

from scripts.bdbv_common import DATA_ROOT, GOLDEN_TEST_OUTBREAK, load_manifest


def _group_stats(split_dir: Path) -> dict:
    if not split_dir.is_dir():
        return {"n_groups": 0, "n_seqs": 0, "outbreaks": {}}
    outbreaks: dict[str, int] = {}
    n_seq = 0
    n_groups = 0
    for gf in sorted(split_dir.glob("*_group_*.fasta")):
        n_groups += 1
        csv_path = gf.with_suffix(".csv")
        if csv_path.is_file():
            import csv

            with open(csv_path, newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    oid = row.get("outbreak_id", "?")
                    outbreaks[oid] = outbreaks.get(oid, 0) + 1
                    n_seq += 1
        else:
            n_seq += sum(1 for _ in SeqIO.parse(gf, "fasta"))
    return {"n_groups": n_groups, "n_seqs": n_seq, "outbreaks": outbreaks}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/filo_l")
    ap.add_argument("--out", default=str(DATA_ROOT / "FILO_SPLIT_AUDIT.json"))
    args = ap.parse_args()

    base = ROOT / args.data
    proto_path = base / "SPLIT_PROTOCOL.json"
    proto = json.loads(proto_path.read_text()) if proto_path.is_file() else {}

    manifest = load_manifest()
    with_outbreak = sum(1 for m in manifest if m.get("outbreak_id"))
    frac = with_outbreak / max(len(manifest), 1)

    report = {
        "data_dir": str(base),
        "protocol": proto,
        "manifest_n": len(manifest),
        "manifest_outbreak_labeled_frac": round(frac, 3),
        "bands": {},
        "ready_to_train": False,
        "blockers": [],
    }
    for split in ("train", "val", "test"):
        report["bands"][split] = _group_stats(base / split)

    train_g = report["bands"]["train"]["n_groups"]
    val_g = report["bands"]["val"]["n_groups"]
    test_g = report["bands"]["test"]["n_groups"]
    train_seq = report["bands"]["train"]["n_seqs"]

    if frac < 0.5:
        report["blockers"].append("outbreak metadata <50% of manifest — re-run Nextstrain ingest")
    if train_g < 3:
        report["blockers"].append(f"train has only {train_g} groups (want >=3)")
    if val_g == 0:
        report["blockers"].append("val empty — will use test for early stopping")
    if test_g == 0:
        report["blockers"].append("test empty — no BDBV 2026 holdout")
    if train_seq < 100:
        report["blockers"].append(f"train only {train_seq} seqs — broaden NCBI/Nextstrain ingest")

    train_golden = report["bands"]["train"]["outbreaks"].get(GOLDEN_TEST_OUTBREAK, 0)
    val_golden = report["bands"]["val"]["outbreaks"].get(GOLDEN_TEST_OUTBREAK, 0)
    test_golden = report["bands"]["test"]["outbreaks"].get(GOLDEN_TEST_OUTBREAK, 0)
    report["golden_test_outbreak"] = GOLDEN_TEST_OUTBREAK
    report["golden_test_in_train"] = train_golden
    report["golden_test_in_val"] = val_golden
    report["golden_test_in_test"] = test_golden
    if train_golden > 0:
        report["blockers"].append(
            f"LEAK: {train_golden} golden-test ({GOLDEN_TEST_OUTBREAK}) seqs in train"
        )
    if val_golden > 0:
        report["blockers"].append(
            f"LEAK: {val_golden} golden-test ({GOLDEN_TEST_OUTBREAK}) seqs in val"
        )
    if test_golden == 0:
        report["blockers"].append(
            f"test missing golden holdout {GOLDEN_TEST_OUTBREAK} — run pathoplexus ingest"
        )

    report["ready_to_train"] = (
        train_g >= 3
        and test_g >= 1
        and train_seq >= 100
        and train_golden == 0
        and val_golden == 0
        and test_golden > 0
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
