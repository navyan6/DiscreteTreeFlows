#!/usr/bin/env python3
"""Audit BDBV temporal vs geographic splits; recommend primary split."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO

from scripts.bdbv_common import DATA_ROOT, parse_date


def _load_split(base: Path) -> dict:
    proto = base / "SPLIT_PROTOCOL.json"
    if proto.is_file():
        return json.loads(proto.read_text())
    return {}


def _count_groups(split_dir: Path) -> tuple[int, int]:
    if not split_dir.is_dir():
        return 0, 0
    fastas = list(split_dir.glob("*_group_*.fasta"))
    n_seq = 0
    for fp in fastas:
        n_seq += sum(1 for _ in SeqIO.parse(fp, "fasta"))
    return len(fastas), n_seq


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--temporal", default="data/bdbv_temporal")
    ap.add_argument("--geo", default="data/bdbv_geo")
    ap.add_argument("--out", default=str(DATA_ROOT / "SPLIT_AUDIT.json"))
    args = ap.parse_args()

    report = {"temporal": {}, "geo": {}, "recommendation": "temporal"}
    for name, rel in (("temporal", args.temporal), ("geo", args.geo)):
        base = ROOT / rel
        proto = _load_split(base)
        bands = {}
        for split in ("train", "val", "test"):
            ng, ns = _count_groups(base / split)
            bands[split] = {"n_groups": ng, "n_seqs": ns}
        report[name] = {"protocol": proto, "bands": bands}

    t_test = report["temporal"]["bands"].get("test", {}).get("n_seqs", 0)
    g_test = report["geo"]["bands"].get("test", {}).get("n_seqs", 0)
    if t_test == 0 and g_test > 0:
        report["recommendation"] = "geo"
    elif t_test > 0:
        report["recommendation"] = "temporal"
    else:
        report["recommendation"] = "need_more_data"

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
