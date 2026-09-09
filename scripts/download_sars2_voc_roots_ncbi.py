#!/usr/bin/env python3
"""
Pull SARS-CoV-2 genomes for VOC-origin populations that are NOT in the TreeSBM
training corpus, so we can build fresh trees whose ROOTS the model has never seen.

Standard applied here (see results/voc_threat_panel/LEAKAGE.md): the model may
have learned escape trajectories from other populations -- shared mutations are
expected and fine -- but it must never have seen this root or this tree. So every
accession is checked against the 110,122 accessions already placed in a
data/covid/{train,val,test} tree.

Selection is by Pango lineage (NCBI Datasets v2alpha `dataset_report` supports
filter.pangolin_classification reliably; filter.geo_location currently 500s, so
country is filtered client-side off isolate name / location fields).

A usable tree needs both the pre-emergence background and the variant itself, so
each population lists `background` lineages (the root context) alongside the
`target` lineages (the VOC being forecast).

Examples:
  # what would be pulled, no network writes
  python scripts/download_sars2_voc_roots_ncbi.py --population usa_epsilon --dry-run

  # real pull
  python scripts/download_sars2_voc_roots_ncbi.py --population usa_epsilon --per-lineage 900
  python scripts/download_sars2_voc_roots_ncbi.py --population all
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# This network sits behind a TLS-intercepting proxy whose root CA is in the macOS
# keychain but not in certifi's bundle. scripts/make_ca_bundle.sh writes a merged
# bundle; pick it up automatically so plain `python scripts/...` just works.
_CA = ROOT / ".certs" / "system_roots.pem"
if _CA.exists() and "SSL_CERT_FILE" not in os.environ:
    os.environ["SSL_CERT_FILE"] = str(_CA)
    os.environ.setdefault("REQUESTS_CA_BUNDLE", str(_CA))

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Entrez.email = "treesbm@example.com"
Entrez.tool = "TreeSBM_voc_roots"

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/virus/taxon/2697049/dataset_report"
SEEN_DEFAULT = ROOT / "results" / "voc_threat_panel" / "seen_accessions.txt"

# country_match is substring-tested against isolate name + location fields.
POPULATIONS = {
    "usa_epsilon": {
        "country": "USA",
        "country_match": ["/USA/", "USA"],
        "date_min": "2020-05-01",
        "date_max": "2021-06-30",
        "background": ["B.1", "B.1.2", "B.1.243"],
        "target": ["B.1.427", "B.1.429"],
        "voc_id": "Epsilon",
        "note": "North America is entirely absent from the training corpus.",
    },
    "usa_iota": {
        "country": "USA",
        "country_match": ["/USA/", "USA"],
        "date_min": "2020-06-01",
        "date_max": "2021-07-31",
        "background": ["B.1", "B.1.2"],
        "target": ["B.1.526"],
        "voc_id": "Iota",
        "note": "North America is entirely absent from the training corpus.",
    },
    "south_africa_beta": {
        "country": "South Africa",
        "country_match": ["/ZAF/", "South Africa"],
        "date_min": "2020-06-01",
        "date_max": "2021-04-30",
        "background": ["B.1", "B.1.1", "C.1", "B.1.1.54", "B.1.1.56"],
        "target": ["B.1.351"],
        "voc_id": "Beta",
        "note": "Country is in train; these accessions are not. New root only.",
    },
    "india_delta": {
        "country": "India",
        "country_match": ["/IND/", "India"],
        "date_min": "2020-09-01",
        "date_max": "2021-07-31",
        "background": ["B.1", "B.1.1", "B.1.36", "B.1.1.306", "B.1.36.29"],
        "target": ["B.1.617.2", "B.1.617.1"],
        "voc_id": "Delta",
        "note": "Country is in train; these accessions are not. New root only.",
    },
    "uk_alpha": {
        "country": "United Kingdom",
        "country_match": ["/GBR/", "United Kingdom", "/England/"],
        "date_min": "2020-06-01",
        "date_max": "2021-03-31",
        "background": ["B.1", "B.1.1"],
        "target": ["B.1.1.7"],
        "voc_id": "Alpha",
        "note": "Country is in train (206 seqs); these accessions are not.",
    },
    "peru_lambda": {
        "country": "Peru",
        "country_match": ["/PER/", "Peru"],
        "date_min": "2020-06-01",
        "date_max": "2021-08-31",
        "background": ["B.1", "B.1.1", "B.1.111", "C.14", "B.1.1.485"],
        "target": ["C.37"],
        "voc_id": "Lambda",
        "note": "Country is in train (165 seqs); these accessions are not.",
    },
    "colombia_mu": {
        "country": "Colombia",
        "country_match": ["/COL/", "Colombia"],
        "date_min": "2020-06-01",
        "date_max": "2021-09-30",
        "background": ["B.1", "B.1.1", "B.1.111", "B.1.420", "B.1.1.348"],
        "target": ["B.1.621", "B.1.621.1"],
        "voc_id": "Mu",
        "note": "Country is in train (29 seqs); these accessions are not.",
    },
}


def _get(url: str, tries: int = 4) -> dict:
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as fh:
                return json.loads(fh.read().decode())
        except Exception as exc:  # noqa: BLE001 - transient 500s are routine here
            print(f"    retry {attempt + 1}/{tries}: {exc}", flush=True)
            time.sleep(2.0 * (attempt + 1))
    return {}


def _loc_blob(rep: dict) -> str:
    iso = rep.get("isolate") or {}
    loc = rep.get("location") or {}
    return " ".join(
        str(x)
        for x in (
            iso.get("name"),
            loc.get("geographic_location"),
            loc.get("geographic_region"),
        )
        if x
    )


def collect_lineage(
    lineage: str,
    country_match: list[str],
    date_min: str,
    date_max: str,
    cap: int,
    seen: set[str],
    max_pages: int = 60,
    min_length: int = 29000,
    page_size: int = 1000,
) -> list[dict]:
    """Page dataset_report for one Pango lineage, keeping in-country, in-window hits."""
    out: list[dict] = []
    token = None
    pages = 0
    while len(out) < cap and pages < max_pages:
        params = {
            "filter.pangolin_classification": lineage,
            "page_size": str(page_size),
        }
        if token:
            params["page_token"] = token
        data = _get(f"{API}?{urllib.parse.urlencode(params)}")
        reps = data.get("reports") or []
        if not reps:
            break
        pages += 1
        for rep in reps:
            acc = (rep.get("accession") or "").split(".")[0]
            if not acc or acc in seen:
                continue
            iso = rep.get("isolate") or {}
            date = (iso.get("collection_date") or "").strip()
            if len(date) < 7 or not (date_min <= date <= date_max):
                continue
            blob = _loc_blob(rep)
            if not any(m in blob for m in country_match):
                continue
            # NCBI flags most SARS-CoV-2 submissions PARTIAL for missing UTR ends,
            # including 29.5kb genomes, so gate on length instead of completeness.
            if (rep.get("length") or 0) < min_length:
                continue
            out.append(
                {
                    "accession": rep.get("accession"),
                    "acc_base": acc,
                    "date": date,
                    "pango": (rep.get("virus") or {}).get("pangolin_classification", ""),
                    "isolate": iso.get("name", ""),
                    "length": rep.get("length"),
                    "completeness": rep.get("completeness", ""),
                }
            )
            if len(out) >= cap:
                break
        token = data.get("next_page_token")
        if not token:
            break
    return out


def fetch_fasta(accs: list[str], batch: int = 100, sleep: float = 0.34) -> dict[str, str]:
    seqs: dict[str, str] = {}
    for i in range(0, len(accs), batch):
        chunk = accs[i : i + batch]
        for attempt in range(3):
            try:
                handle = Entrez.efetch(
                    db="nucleotide", id=",".join(chunk), rettype="fasta", retmode="text"
                )
                for rec in SeqIO.parse(handle, "fasta"):
                    seqs[rec.id.split(".")[0]] = str(rec.seq)
                handle.close()
                break
            except Exception as exc:  # noqa: BLE001
                print(f"    efetch retry {attempt + 1}: {exc}", flush=True)
                time.sleep(2.0 * (attempt + 1))
        print(f"    fetched {len(seqs)}/{len(accs)}", flush=True)
        time.sleep(sleep)
    return seqs


def run_population(key: str, args: argparse.Namespace, seen: set[str]) -> dict:
    pop = POPULATIONS[key]
    print(f"\n=== {key}  ({pop['country']}, {pop['date_min']}..{pop['date_max']}) ===")
    records: list[dict] = []
    for role in ("background", "target"):
        cap = args.per_lineage if role == "background" else args.per_lineage_target
        for lin in pop[role]:
            hits = collect_lineage(
                lin, pop["country_match"], pop["date_min"], pop["date_max"], cap, seen,
                max_pages=args.max_pages, min_length=args.min_length,
            )
            for h in hits:
                h["role"] = role
            records.extend(hits)
            print(f"  {role:<10} {lin:<12} -> {len(hits)} unseen in-country hits")

    uniq: dict[str, dict] = {}
    for r in records:
        uniq.setdefault(r["acc_base"], r)
    records = sorted(uniq.values(), key=lambda r: r["date"])
    n_t = sum(1 for r in records if r["role"] == "target")
    print(f"  total {len(records)} unique ({n_t} target-lineage, {len(records) - n_t} background)")

    out_dir = ROOT / args.out_base / "raw"
    out_dir.mkdir(parents=True, exist_ok=True)
    meta_path = out_dir / f"{key}_metadata.csv"
    with meta_path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["accession", "acc_base", "date", "pango", "role", "isolate",
                        "length", "completeness"],
        )
        w.writeheader()
        w.writerows(records)
    print(f"  wrote {meta_path.relative_to(ROOT)}")

    if args.dry_run or not records:
        return {"population": key, "n_records": len(records), "n_target": n_t, "fasta": None}

    seqs = fetch_fasta([r["accession"] for r in records])
    fasta_path = out_dir / f"{key}_covid_seqs.fasta"
    written = []
    for r in records:
        seq = seqs.get(r["acc_base"])
        if not seq:
            continue
        # Pipeline header contract (prepare_covid_geo.parse_header):
        # split('|') -> [acc, desc, date, length, country]
        desc = f"|{r['pango']} {r['isolate']}|{r['date']}|{len(seq)}|{pop['country']}"
        written.append(SeqRecord(seq=Seq(seq), id=r["acc_base"], description=desc))
    SeqIO.write(written, fasta_path, "fasta")
    print(f"  wrote {len(written)} seqs -> {fasta_path.relative_to(ROOT)}")
    return {
        "population": key,
        "voc_id": pop["voc_id"],
        "country": pop["country"],
        "n_records": len(records),
        "n_target": n_t,
        "n_written": len(written),
        "fasta": str(fasta_path.relative_to(ROOT)),
        "metadata": str(meta_path.relative_to(ROOT)),
        "note": pop["note"],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--population", default="usa_epsilon", help="key from POPULATIONS, or 'all'")
    ap.add_argument("--per-lineage", type=int, default=700, help="cap per background lineage")
    ap.add_argument("--per-lineage-target", type=int, default=700, help="cap per target lineage")
    ap.add_argument("--max-pages", type=int, default=60,
                    help="pages of 1000 per lineage; raise for sparse country/lineage pairs")
    ap.add_argument("--min-length", type=int, default=29000,
                    help="min genome length; NCBI marks 29.5kb genomes PARTIAL")
    ap.add_argument("--out-base", default="data/covid_voc_roots")
    ap.add_argument("--seen-accessions", type=Path, default=SEEN_DEFAULT)
    ap.add_argument("--dry-run", action="store_true", help="screen only, no sequence download")
    args = ap.parse_args()

    seen: set[str] = set()
    if args.seen_accessions.exists():
        seen = {
            line.strip().split(".")[0]
            for line in args.seen_accessions.read_text().splitlines()
            if line.strip()
        }
    print(f"excluding {len(seen)} accessions already used in a data/covid tree")

    keys = list(POPULATIONS) if args.population == "all" else [args.population]
    unknown = [k for k in keys if k not in POPULATIONS]
    if unknown:
        ap.error(f"unknown population(s) {unknown}; choose from {list(POPULATIONS)}")

    summary = [run_population(k, args, seen) for k in keys]
    out = ROOT / args.out_base / "raw" / "pull_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)}")
    print("Next: scripts/prepare_covid_voc_roots.py, then the Betty spike+tree pipeline.")


if __name__ == "__main__":
    main()
