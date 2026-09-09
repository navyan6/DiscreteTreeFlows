#!/usr/bin/env python3
"""
Download complete filovirus genomes from NCBI nucleotide database.

Queries per pathogenic Orthoebolavirus species, fetches GenBank records,
writes per-species FASTA + unified manifest.json + DATA_AUDIT.md.

Usage:
  python scripts/download_ebolavirus_ncbi.py
  python scripts/download_ebolavirus_ncbi.py --max-per-species 500 --retmax 5000
  python scripts/download_ebolavirus_ncbi.py --dry-run
"""

from __future__ import annotations

import argparse
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import Entrez, SeqIO

from scripts.bdbv_common import (
    DATA_ROOT,
    SPECIES,
    dedupe_manifest,
    enrich_manifest_outbreaks,
    infer_outbreak_id,
    load_manifest,
    normalize_country,
    parse_date,
    save_manifest,
    write_data_audit,
)

Entrez.email = "treesbm@example.com"


def _parse_genbank_date(raw: str) -> str:
    raw = (raw or "").strip()
    for fmt in ("%d-%b-%Y", "%Y-%m-%d", "%Y/%m/%d", "%Y"):
        try:
            from datetime import datetime

            dt = datetime.strptime(raw[: len(fmt.replace("%", "0"))], fmt)
            return dt.strftime("%Y-%m-%d")
        except ValueError:
            continue
    m = re.search(r"(\d{4})", raw)
    return f"{m.group(1)}-01-01" if m else ""


def _country_from_record(rec) -> str:
    for key in ("geo_loc_name", "country", "Isolation-source"):
        if key in rec.annotations:
            return normalize_country(str(rec.annotations[key]))
    desc = rec.description or ""
    for part in desc.split("|"):
        if ":" in part and len(part) < 80:
            return normalize_country(part.split(":")[0])
    return "unknown"


def _host_from_record(rec) -> str:
    return str(rec.annotations.get("host", rec.annotations.get("organism", "")))


def _queries_for_species(sp_key: str, org: str) -> list[str]:
    if sp_key == "marv":
        return [
            "Marburgvirus[Organism] AND 17000:25000[SLEN]",
            "txid11269[Organism:exp] AND 17000:25000[SLEN]",
            '"Marburg marburgvirus"[Organism] AND 17000:25000[SLEN]',
        ]
    return [
        (
            f'"{org}"[Organism] AND 17000:25000[SLEN] '
            f'AND ("complete"[Title] OR "genome"[Title] OR "complete genome"[Title])'
        )
    ]


def _esearch_ids_union(queries: list[str], retmax: int) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for query in queries:
        print(f"  query: {query}")
        for uid in esearch_ids(query, retmax):
            if uid not in seen:
                seen.add(uid)
                ordered.append(uid)
    return ordered


def esearch_ids(query: str, retmax: int) -> list[str]:
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])


def efetch_genbank(ids: list[str]) -> list:
    if not ids:
        return []
    records = []
    batch = 200
    for i in range(0, len(ids), batch):
        chunk = ids[i : i + batch]
        handle = Entrez.efetch(
            db="nucleotide", id=chunk, rettype="gb", retmode="text"
        )
        records.extend(SeqIO.parse(handle, "genbank"))
        handle.close()
        time.sleep(0.34)
    return records


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--retmax", type=int, default=8000, help="Max IDs per species query")
    ap.add_argument("--max-per-species", type=int, default=0, help="Cap fetched records (0=all)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--email", default="treesbm@example.com")
    ap.add_argument(
        "--species",
        action="append",
        default=[],
        help="Only fetch these species keys (repeatable). Default: all in SPECIES.",
    )
    args = ap.parse_args()
    Entrez.email = args.email

    species_keys = args.species or list(SPECIES.keys())

    out_dir = DATA_ROOT / "raw" / "ncbi"
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_dir = DATA_ROOT / "references"
    ref_dir.mkdir(parents=True, exist_ok=True)

    manifest_new: list[dict] = []
    for sp_key, meta in SPECIES.items():
        if sp_key not in species_keys:
            continue
        org = meta["organism"]
        queries = _queries_for_species(sp_key, org)
        print(f"[{sp_key}] {len(queries)} query/queries")
        ids = _esearch_ids_union(queries, args.retmax)
        print(f"  esearch union: {len(ids)} ids")
        if args.dry_run:
            continue
        records = efetch_genbank(ids)
        if args.max_per_species > 0:
            records = records[: args.max_per_species]
        sp_fasta = out_dir / f"{sp_key}_genomes.fasta"
        n_written = 0
        with open(sp_fasta, "w") as out_f:
            for rec in records:
                if len(rec.seq) < 15000:
                    continue
                acc = rec.id.split(".")[0]
                date_raw = _parse_genbank_date(
                    rec.annotations.get("collection_date", "")
                    or rec.annotations.get("date", "")
                )
                augur_date, _, year = parse_date(date_raw)
                country = _country_from_record(rec)
                host = _host_from_record(rec)
                out_f.write(f">{acc},{augur_date},{country},{sp_key},ncbi\n")
                out_f.write(f"{rec.seq}\n")
                n_written += 1
                manifest_new.append(
                    {
                        "accession": acc,
                        "species": sp_key,
                        "date": augur_date,
                        "year": year,
                        "country": country,
                        "source": "ncbi",
                        "length": len(rec.seq),
                        "host": host,
                        "fasta": str(sp_fasta.relative_to(ROOT)),
                        "outbreak_id": infer_outbreak_id(
                            {
                                "accession": acc,
                                "species": sp_key,
                                "year": year,
                                "country": country,
                                "source": "ncbi",
                            }
                        ),
                    }
                )
        print(f"  wrote {n_written} genomes -> {sp_fasta}")

    if args.dry_run:
        print("Dry run — no files written.")
        return

    # Fetch reference genomes
    for acc in ("NC_014373.1", "FJ217161.1", "KM034562.1", "NC_001608.1"):
        try:
            recs = efetch_genbank(esearch_ids(acc, 1))
            if recs:
                SeqIO.write(recs[0], ref_dir / f"{recs[0].id}.gb", "genbank")
                SeqIO.write(recs[0], ref_dir / f"{recs[0].id}.fasta", "fasta")
                print(f"Reference: {ref_dir / recs[0].id}.fasta")
        except Exception as e:
            print(f"WARNING: reference {acc}: {e}")

    existing = load_manifest()
    if args.species:
        drop = set(species_keys)
        existing = [
            r
            for r in existing
            if not (r.get("source") == "ncbi" and r.get("species") in drop)
        ]
    else:
        existing = [r for r in existing if r.get("source") != "ncbi"]
    merged = dedupe_manifest(existing + manifest_new)
    save_manifest(merged)
    write_data_audit(merged)
    print(f"Manifest: {len(merged)} records -> {DATA_ROOT / 'manifest.json'}")


if __name__ == "__main__":
    main()
