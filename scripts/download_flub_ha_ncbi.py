#!/usr/bin/env python3
"""
Download human Influenza B HA (segment 4) from NCBI for pan-flu training.

Pulls Victoria + Yamagata lineages (post-2009 pandemic window), writes:
  data/flub/train/flubtrain.fasta
  data/flub/SPLIT_PROTOCOL.json  (inventory only)

Headers: >ACCESSION,DATE,lineage={victoria|yamagata|unknown},subtype=flub

Usage:
  python scripts/download_flub_ha_ncbi.py
  python scripts/download_flub_ha_ncbi.py --dry-run --max-records 500
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Entrez.email = "treesbm@example.com"

QUERIES = [
    (
        "victoria",
        '"Influenza B virus"[Organism] AND (Victoria[All Fields] OR "B/Victoria"[All Fields]) '
        "AND 1400:2000[SLEN] AND 2009:2026[dp]",
    ),
    (
        "yamagata",
        '"Influenza B virus"[Organism] AND (Yamagata[All Fields] OR "B/Yamagata"[All Fields]) '
        "AND 1400:2000[SLEN] AND 2009:2020[dp]",
    ),
    (
        "unknown",
        '"Influenza B virus"[Organism] AND HA[Title] AND 1400:2000[SLEN] AND 2009:2026[dp]',
    ),
]


def _parse_date(raw: str) -> str:
    raw = (raw or "").strip()
    for fmt in ("%d-%b-%Y", "%Y-%m-%d", "%Y/%m/%d", "%Y"):
        try:
            from datetime import datetime

            dt = datetime.strptime(raw[: len(fmt.replace("%", "0"))], fmt)
            return dt.strftime("%Y-%m-%d")
        except ValueError:
            continue
    m = re.search(r"(\d{4})", raw)
    return f"{m.group(1)}-07-01" if m else "2000-01-01"


def _infer_lineage(rec, default: str) -> str:
    blob = " ".join(
        [
            str(rec.description or ""),
            str(rec.annotations.get("organism", "")),
            str(rec.annotations.get("source", "")),
        ]
    ).lower()
    if "yamagata" in blob or "/b/yamagata" in blob:
        return "yamagata"
    if "victoria" in blob or "/b/victoria" in blob:
        return "victoria"
    return default if default != "unknown" else "unknown"


def _is_ha_segment(rec) -> bool:
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        prod = feat.qualifiers.get("product", [""])[0].lower()
        note = " ".join(feat.qualifiers.get("note", [])).lower()
        if "hemagglutinin" in prod or "ha protein" in prod or "segment 4" in note:
            return True
    desc = (rec.description or "").lower()
    return "hemagglutinin" in desc or " segment 4 " in desc or " ha " in desc


def _extract_ha_nt(rec) -> str | None:
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        prod = feat.qualifiers.get("product", [""])[0].lower()
        if "hemagglutinin" in prod or "ha" == prod.strip():
            return str(feat.extract(rec.seq))
    if 1400 <= len(rec.seq) <= 2000:
        return str(rec.seq)
    return None


def esearch_ids(query: str, retmax: int) -> list[str]:
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])


def efetch_genbank(ids: list[str]) -> list:
    out = []
    for i in range(0, len(ids), 200):
        chunk = ids[i : i + 200]
        handle = Entrez.efetch(db="nucleotide", id=chunk, rettype="gb", retmode="text")
        out.extend(SeqIO.parse(handle, "genbank"))
        handle.close()
        time.sleep(0.34)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-base", default="data/flub")
    ap.add_argument("--retmax", type=int, default=8000)
    ap.add_argument("--max-records", type=int, default=12000)
    ap.add_argument("--min-year", type=int, default=2009)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    seen: set[str] = set()
    rows: list[tuple[str, str, str, str, str]] = []
    lineage_counts: Counter = Counter()

    for default_lineage, query in QUERIES:
        print(f"=== {default_lineage}: {query}")
        ids = esearch_ids(query, args.retmax)
        print(f"  ids={len(ids)}")
        if args.dry_run:
            continue
        for rec in efetch_genbank(ids):
            acc = rec.id.split(".")[0]
            if acc in seen:
                continue
            if not _is_ha_segment(rec):
                continue
            nt = _extract_ha_nt(rec)
            if not nt or "N" in nt.upper()[:30]:
                continue
            try:
                aa = str(Seq(nt).translate(to_stop=True))
            except Exception:
                continue
            if len(aa) < 500:
                continue
            date = _parse_date(rec.annotations.get("date", ""))
            try:
                if int(date[:4]) < args.min_year:
                    continue
            except ValueError:
                pass
            lineage = _infer_lineage(rec, default_lineage)
            seen.add(acc)
            rows.append((acc, date, aa[:566], lineage, "flub"))
            lineage_counts[lineage] += 1
            if len(rows) >= args.max_records:
                break
        if len(rows) >= args.max_records:
            break

    if args.dry_run:
        print("Dry run — no FASTA written")
        return

    base = ROOT / args.out_base
    train_dir = base / "train"
    train_dir.mkdir(parents=True, exist_ok=True)
    out_fa = train_dir / "flubtrain.fasta"
    with open(out_fa, "w") as ff:
        for acc, date, aa, lineage, subtype in rows:
            ff.write(f">{acc},{date},lineage={lineage},subtype={subtype}\n{aa}\n")

    protocol = {
        "dataset": "flub",
        "split_type": "inventory",
        "source": "NCBI nucleotide (segment 4 / HA CDS)",
        "lineages": dict(lineage_counts),
        "n_seqs": len(rows),
        "min_year": args.min_year,
        "out_fasta": str(out_fa.relative_to(ROOT)),
        "note": "Use prepare_panflu_pool.py --include-flub after ingest.",
    }
    (base / "SPLIT_PROTOCOL.json").write_text(json.dumps(protocol, indent=2) + "\n")
    print(f"Wrote {out_fa} n={len(rows)} lineages={dict(lineage_counts)}")


if __name__ == "__main__":
    main()
