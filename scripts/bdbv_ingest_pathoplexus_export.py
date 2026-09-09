#!/usr/bin/env python3
"""
Ingest Pathoplexus BDBV export (metadata TSV + aligned-nuc FASTA).

Accepts plain or .gz files. Writes:
  data/bdbv/raw/pathoplexus_2026/bdbv_genomes.fasta
  merges manifest (source=pathoplexus, outbreak_id bdbv_2026 for 2026 band)

Usage:
  python scripts/bdbv_ingest_pathoplexus_export.py
  python scripts/bdbv_ingest_pathoplexus_export.py \\
    --metadata data/bdbv/ebola-bdbv_metadata_2026-08-31T0254.tsv.gz \\
    --fasta data/bdbv/ebola-bdbv_aligned-nuc_2026-08-31T0254.fasta.gz
  python scripts/bdbv_ingest_pathoplexus_export.py --resplit-only
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scripts.bdbv_common import (
    DATA_ROOT,
    GOLDEN_TEST_OUTBREAK,
    dedupe_manifest,
    enrich_manifest_outbreaks,
    format_fasta_header,
    infer_outbreak_id,
    load_manifest,
    normalize_country,
    parse_date,
    save_manifest,
    slug_outbreak_id,
    write_data_audit,
)

INSDC_RE = re.compile(r"^[A-Z]{1,2}\d+[A-Z]?\d*$")
DEFAULT_META_GLOB = "ebola-bdbv_metadata_*.tsv*"
DEFAULT_FASTA_GLOB = "ebola-bdbv_aligned-nuc_*.fasta*"
OUT_DIR = DATA_ROOT / "raw" / "pathoplexus_2026"
GENOME_FASTA = OUT_DIR / "bdbv_genomes.fasta"


def _open_text(path: Path):
    if path.suffix == ".gz" or path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, newline="", encoding="utf-8")


def _find_export_file(base: Path, pattern: str) -> Path | None:
    matches = sorted(base.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0] if matches else None


def _outbreak_id_from_label(outbreak: str, year: int | None) -> str:
    ob = (outbreak or "").strip()
    if ob.lower() == "bdbv-2026":
        return "bdbv_2026"
    if ob.lower() == "bdbv-2012":
        return "bdbv_2012"
    if ob.lower() == "bdbv-2007":
        return "bdbv_2007"
    if ob and ob.lower() not in ("?", "unknown", "unassigned"):
        return slug_outbreak_id(f"bdbv_{ob}")
    return infer_outbreak_id({"species": "bdbv", "year": year, "outbreak": ob})


def _canonical_accession(pp_id: str, insdc: str) -> str:
    insdc = (insdc or "").strip().split(".")[0]
    if insdc and INSDC_RE.match(insdc):
        return insdc
    return pp_id.split(".")[0]


def load_metadata(path: Path) -> dict[str, dict]:
    meta_by_ppx: dict[str, dict] = {}
    with _open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ppx = (row.get("accessionVersion") or "").strip()
            if not ppx:
                continue
            ppx_key = ppx.split(".")[0]
            date_raw = (row.get("sampleCollectionDate") or "").strip()
            augur_date, _, year = parse_date(date_raw)
            country = normalize_country(row.get("geoLocCountry") or row.get("geoLocAdmin1") or "")
            outbreak = (row.get("outbreak") or "").strip()
            insdc = (row.get("insdcAccessionFull") or "").strip()
            acc = _canonical_accession(ppx, insdc)
            outbreak_id = _outbreak_id_from_label(outbreak, year)
            # Golden eval band only — historical BDBV stays on NCBI/Nextstrain train paths.
            if outbreak_id != GOLDEN_TEST_OUTBREAK:
                continue
            meta_by_ppx[ppx_key] = {
                "accession": acc,
                "pathoplexus_id": ppx_key,
                "species": "bdbv",
                "date": augur_date,
                "year": year,
                "country": country,
                "source": "pathoplexus",
                "host": row.get("hostNameScientific") or "Homo sapiens",
                "outbreak": outbreak,
                "outbreak_id": outbreak_id,
                "golden_test": True,
                "division": row.get("geoLocAdmin1") or "",
                "data_use_terms": row.get("dataUseTerms") or "",
                "insdc_accession": insdc.split(".")[0] if insdc else "",
            }
    return meta_by_ppx


def _clean_genome(seq: str) -> str:
    return re.sub(r"[^ACGTN]", "", seq.upper())


def ingest_export(metadata: Path, fasta: Path) -> list[dict]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta_by_ppx = load_metadata(metadata)
    print(f"Metadata: {len(meta_by_ppx)} Pathoplexus rows from {metadata.name}")

    records: list[dict] = []
    n_missing_meta = 0
    with open(GENOME_FASTA, "w") as out_f:
        for rec in SeqIO.parse(_open_text(fasta), "fasta"):
            ppx_key = rec.id.split(".")[0].split()[0]
            m = meta_by_ppx.get(ppx_key)
            if m is None:
                n_missing_meta += 1
                continue
            seq = _clean_genome(str(rec.seq))
            if len(seq) < 15000:
                continue
            acc = m["accession"]
            hdr = format_fasta_header(
                acc, m["date"], m["country"], m["species"], m["source"]
            )
            out_f.write(f"{hdr}\n{seq}\n")
            row = dict(m)
            row["length"] = len(seq)
            row["fasta"] = str(GENOME_FASTA.relative_to(ROOT))
            records.append(row)

    print(
        f"Wrote {len(records)} genomes -> {GENOME_FASTA} "
        f"(missing meta for {n_missing_meta} FASTA records)"
    )
    by_outbreak: dict[str, int] = {}
    for r in records:
        by_outbreak[r["outbreak_id"]] = by_outbreak.get(r["outbreak_id"], 0) + 1
    print("By outbreak_id:", dict(sorted(by_outbreak.items(), key=lambda x: -x[1])))
    return records


def merge_manifest(new_records: list[dict]) -> tuple[int, int]:
    existing = load_manifest()
    before = len(existing)
    # Drop prior pathoplexus rows so re-ingest is idempotent
    existing = [r for r in existing if r.get("source") != "pathoplexus"]
    merged = dedupe_manifest(existing + new_records)
    merged = enrich_manifest_outbreaks(merged)
    save_manifest(merged)
    write_data_audit(merged)
    audit = OUT_DIR / "pathoplexus_ingest_audit.json"
    audit.write_text(
        json.dumps(
            {
                "pathoplexus_genomes": len(new_records),
                "manifest_total": len(merged),
                "test_outbreak_bdbv_2026": sum(
                    1 for r in new_records if r.get("outbreak_id") == "bdbv_2026"
                ),
            },
            indent=2,
        )
        + "\n"
    )
    return before, len(merged)


def resplit_filo(out_base: str = "data/filo_l", group_size: int = 80, min_group: int = 5) -> None:
    import subprocess

    py = sys.executable
    steps = [
        [py, "scripts/bdbv_extract_l.py", "--pathoplexus-only"],
        [py, "scripts/bdbv_define_l_window.py"],
        [
            py,
            "scripts/prepare_filo_outbreak.py",
            "--out-base",
            out_base,
            "--group-size",
            str(group_size),
            "--min-group",
            str(min_group),
        ],
        [py, "scripts/audit_filo_splits.py", "--data", out_base],
    ]
    for cmd in steps:
        print(f"\n>>> {' '.join(cmd)}", flush=True)
        subprocess.run(cmd, cwd=ROOT, check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", default="", help="Pathoplexus metadata TSV (.gz ok)")
    ap.add_argument("--fasta", default="", help="Aligned nucleotide FASTA (.gz ok)")
    ap.add_argument(
        "--search-dir",
        default=str(DATA_ROOT),
        help="Directory to search for export files if --metadata/--fasta omitted",
    )
    ap.add_argument("--resplit-only", action="store_true", help="Skip ingest; re-extract L + filo splits")
    ap.add_argument("--no-resplit", action="store_true", help="Ingest only; do not rebuild filo_l splits")
    ap.add_argument("--out-base", default="data/filo_l")
    ap.add_argument("--group-size", type=int, default=80)
    ap.add_argument("--min-group", type=int, default=5)
    args = ap.parse_args()

    if not args.resplit_only:
        search = Path(args.search_dir)
        meta_path = Path(args.metadata) if args.metadata else _find_export_file(search, DEFAULT_META_GLOB)
        fasta_path = Path(args.fasta) if args.fasta else _find_export_file(search, DEFAULT_FASTA_GLOB)
        if meta_path is None or not meta_path.is_file():
            raise SystemExit(f"Missing Pathoplexus metadata under {search} ({DEFAULT_META_GLOB})")
        if fasta_path is None or not fasta_path.is_file():
            raise SystemExit(f"Missing Pathoplexus FASTA under {search} ({DEFAULT_FASTA_GLOB})")

        records = ingest_export(meta_path, fasta_path)
        before, after = merge_manifest(records)
        print(f"Manifest: {before} -> {after} (+{after - before} net after dedupe)")

    if not args.no_resplit:
        resplit_filo(args.out_base, args.group_size, args.min_group)
        test_dir = ROOT / args.out_base / "test"
        ng = len(list(test_dir.glob("*_group_*.fasta"))) if test_dir.is_dir() else 0
        print(f"\nTest trees ready: {ng} groups under {test_dir}")


if __name__ == "__main__":
    main()
