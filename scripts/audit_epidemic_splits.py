#!/usr/bin/env python3
"""Audit epidemic split dirs (filo / COVID / flu)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.epidemic_split_common import count_groups


def audit_one(data_dir: Path) -> dict:
    proto_path = data_dir / "SPLIT_PROTOCOL.json"
    proto = json.loads(proto_path.read_text()) if proto_path.is_file() else {}
    bands = {}
    for split in ("train", "val", "test"):
        bands[split] = {"dir": str(data_dir / split), **dict(zip(
            ("n_groups", "n_seqs"), count_groups(data_dir / split)
        ))}
    blockers = []
    if bands["train"]["n_groups"] < 2:
        blockers.append("train has <2 groups")
    if bands["val"]["n_groups"] == 0:
        blockers.append("val empty")
    if bands["test"]["n_groups"] == 0:
        blockers.append("test empty")
    return {
        "data_dir": str(data_dir),
        "protocol": proto,
        "bands": bands,
        "blockers": blockers,
        "ready_to_train": len(blockers) == 0 or (
            bands["train"]["n_groups"] >= 2 and bands["test"]["n_groups"] >= 1
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--data",
        action="append",
        default=[],
        help="Data dir (repeatable). Defaults: filo_l, covid_epidemic, h3n2_epidemic, h1n1_epidemic",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    dirs = args.data or [
        "data/filo_l",
        "data/covid_epidemic",
        "data/h3n2_epidemic",
        "data/h1n1_epidemic",
    ]
    report = {d: audit_one(ROOT / d) for d in dirs}
    text = json.dumps(report, indent=2) + "\n"
    if args.out:
        Path(args.out).write_text(text)
    print(text)


if __name__ == "__main__":
    main()
