#!/usr/bin/env python3
"""
Assemble VOC threat case-study folders from existing COVID trees.

Reads candidates JSON (from screen_voc_threat_trees.py) or an explicit
--cases list and copies anc_aa + rooted.nwk (+ meta csv) into:

  results/voc_threat_panel/cases/{voc_id}_g{gid:03d}_{split}/

Also writes cases_manifest.json for betty_submit_voc_threat_panel.sh.

Examples:
  python scripts/prepare_voc_origin_cases.py \\
    --candidates results/voc_threat_panel/candidates_test.json \\
    --candidates results/voc_threat_panel/candidates_train.json \\
    --max-per-voc 2

  python scripts/prepare_voc_origin_cases.py \\
    --case Gamma:test:3 --case Beta:train:302 --case Delta:train:191
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from voc_threat_lib import group_paths, load_panel  # noqa: E402


def parse_case(s: str) -> tuple[str, str, int]:
    # VOC:split:gid
    parts = s.split(":")
    if len(parts) != 3:
        raise SystemExit(f"Bad --case {s!r}; expected VOC:split:gid")
    return parts[0], parts[1], int(parts[2])


def data_dir_for(split: str) -> Path:
    return ROOT / "data" / "covid" / split


def copy_case(voc_id: str, split: str, gid: int, out_root: Path, meta: dict | None) -> dict:
    src = data_dir_for(split)
    paths = group_paths(src, gid)
    if not paths["anc"].exists():
        raise FileNotFoundError(paths["anc"])
    dest = out_root / f"{voc_id}_g{gid:03d}_{split}"
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(paths["anc"], dest / "observed_anc_aa.fasta")
    if paths["nwk"].exists():
        shutil.copy2(paths["nwk"], dest / "observed.nwk")
    for key in ("csv", "csv_train", "csv_val"):
        if paths[key].exists():
            shutil.copy2(paths[key], dest / "meta.csv")
            break
    info = {
        "voc_id": voc_id,
        "split": split,
        "gid": gid,
        "case_dir": str(dest.relative_to(ROOT)),
        "source_data_dir": str(src),
        "n_leaves_hint": meta.get("n_leaves") if meta else None,
        "date_min": meta.get("date_min") if meta else None,
        "date_max": meta.get("date_max") if meta else None,
        "screen": meta,
    }
    (dest / "case.json").write_text(json.dumps(info, indent=2) + "\n")
    return info


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=Path, default=None)
    ap.add_argument("--candidates", type=Path, action="append", default=[])
    ap.add_argument("--case", action="append", default=[], help="VOC:split:gid")
    ap.add_argument("--max-per-voc", type=int, default=2)
    ap.add_argument("--prefer-vocs", type=str, default="")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "voc_threat_panel" / "cases",
    )
    ap.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "results" / "voc_threat_panel" / "cases_manifest.json",
    )
    args = ap.parse_args()

    panel = load_panel(args.panel)
    prefer = [
        x.strip()
        for x in (args.prefer_vocs or ",".join(panel.get("primary_panel", []))).split(",")
        if x.strip()
    ]

    selected: list[tuple[str, str, int, dict | None]] = []
    seen = set()

    for c in args.case:
        voc, split, gid = parse_case(c)
        key = (voc, split, gid)
        if key not in seen:
            selected.append((voc, split, gid, None))
            seen.add(key)

    # From candidate JSONs: take top per VOC
    pooled = []
    for cp in args.candidates:
        blob = json.loads(cp.read_text())
        pooled.extend(blob.get("candidates") or [])
    pooled.sort(key=lambda x: x.get("score", 0), reverse=True)

    per_voc_count: dict[str, int] = {}
    for c in pooled:
        voc = c["voc_id"]
        if prefer and voc not in prefer:
            continue
        n = per_voc_count.get(voc, 0)
        if n >= args.max_per_voc:
            continue
        key = (voc, c["split"], int(c["gid"]))
        if key in seen:
            continue
        selected.append((voc, c["split"], int(c["gid"]), c))
        seen.add(key)
        per_voc_count[voc] = n + 1

    # Seed defaults from panel if still empty for primary VOCs
    if not selected:
        seeds = [
            ("Gamma", "test", 3),
            ("Beta", "train", 302),
            ("Delta", "train", 191),
            ("Omicron_BA1", "train", 309),
            ("Alpha", "train", 323),
            ("Lambda", "train", 289),
            ("Mu", "train", 34),
        ]
        for voc, split, gid in seeds:
            selected.append((voc, split, gid, None))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    cases = []
    for voc, split, gid, meta in selected:
        try:
            info = copy_case(voc, split, gid, args.out_dir, meta)
            cases.append(info)
            print(f"ok {voc} {split} g{gid:03d} -> {info['case_dir']}")
        except FileNotFoundError as e:
            print(f"SKIP {voc} {split} g{gid:03d}: {e}")

    manifest = {
        "n_cases": len(cases),
        "cases": cases,
        "generation_defaults": panel.get("generation_defaults", {}),
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote manifest {args.manifest} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
