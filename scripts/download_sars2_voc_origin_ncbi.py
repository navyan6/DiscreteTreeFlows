#!/usr/bin/env python3
"""
Optional NCBI top-up for sparse VOC-origin Spike AA sets (Mu/Colombia, Lambda/Peru).

Pulls SARS-CoV-2 complete genomes from NCBI Virus/Nucleotide for a country +
date window, writes raw genomes in the same header format as existing
*_covid_seqs.fasta so covid_extract_spike.py / prepare_covid_geo.py can ingest.

NOTE: Spike AA extraction still needs nextclade (Betty). This script only
fetches genomes. Prefer existing train groups when available.

Examples:
  python scripts/download_sars2_voc_origin_ncbi.py --country Colombia \\
    --date-min 2021-01-01 --date-max 2021-08-01 --max-records 500 --dry-run

  python scripts/download_sars2_voc_origin_ncbi.py --country Peru \\
    --date-min 2020-12-01 --date-max 2021-07-01 --max-records 800
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = "treesbm@example.com"
Entrez.tool = "TreeSBM_voc_threat"


COUNTRY_ALIASES = {
    "Colombia": "Colombia",
    "Peru": "Peru",
    "South Africa": "South Africa",
    "United Kingdom": "United Kingdom",
    "UK": "United Kingdom",
    "India": "India",
    "Brazil": "Brazil",
}


def esearch_ids(query: str, retmax: int) -> list[str]:
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y")
    rec = Entrez.read(handle)
    handle.close()
    return list(rec.get("IdList") or [])


def efetch_fasta(ids: list[str], sleep: float = 0.35) -> list:
    out = []
    for i in range(0, len(ids), 50):
        chunk = ids[i : i + 50]
        for attempt in range(3):
            try:
                handle = Entrez.efetch(
                    db="nucleotide", id=",".join(chunk), rettype="fasta", retmode="text"
                )
                out.extend(list(SeqIO.parse(handle, "fasta")))
                handle.close()
                break
            except Exception as e:
                print(f"efetch retry {attempt+1}: {e}", flush=True)
                time.sleep(1.5 * (attempt + 1))
        time.sleep(sleep)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--country", required=True)
    ap.add_argument("--date-min", default="2020-01-01")
    ap.add_argument("--date-max", default="2022-12-31")
    ap.add_argument("--max-records", type=int, default=500)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Default: data/covid/train/ncbi_{country}_covid_seqs.fasta",
    )
    args = ap.parse_args()

    country = COUNTRY_ALIASES.get(args.country, args.country)
    # NCBI Virus-style query for SARS-CoV-2 complete genomes
    query = (
        f'Severe acute respiratory syndrome coronavirus 2[Organism] AND '
        f'complete genome[Title] AND {country}[Country] AND '
        f'{args.date_min}:{args.date_max}[PDAT]'
    )
    print("query:", query, flush=True)
    ids = esearch_ids(query, args.max_records)
    print(f"esearch hits (capped): {len(ids)}", flush=True)
    if args.dry_run or not ids:
        return

    recs = efetch_fasta(ids)
    out = args.out or (
        ROOT / "data" / "covid" / "train" / f"ncbi_{country.replace(' ', '_').lower()}_covid_seqs.fasta"
    )
    out.parent.mkdir(parents=True, exist_ok=True)

    written = []
    for rec in recs:
        acc = rec.id.split(".")[0] + ("." + rec.id.split(".")[1] if "." in rec.id else "")
        # Normalize to pipeline header: ACCESSION |desc|date|length|country
        # Date often absent in fasta-only fetch — leave placeholder for extract step
        desc = (
            f"{rec.id} |{rec.description}|{args.date_min}|{len(rec.seq)}|{country}"
        )
        written.append(SeqRecord(rec.seq, id=rec.id, description=desc))

    SeqIO.write(written, out, "fasta")
    print(f"wrote {len(written)} records -> {out}")
    print("Next on Betty: run covid_extract_spike.py on this file, then rebuild groups.")


if __name__ == "__main__":
    main()
