#!/usr/bin/env python3
"""
Download human influenza PB2 (and optionally NP) CDS from NCBI Influenza DB / GenBank.

Writes AA FASTA under data/<gene_id>/raw/ for later epidemic/season grouping.

Usage:
  python scripts/download_flu_segment_ncbi.py --subtype H3N2 --segment PB2 --retmax 5000
  python scripts/download_flu_segment_ncbi.py --subtype H3N2 --segment NP --retmax 5000
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parent.parent
Entrez.email = "nnori@upenn.edu"

SEG_LEN = {
    "PB2": (2200, 2400),
    "PB1": (2200, 2400),
    "PA": (2100, 2300),
    "NP": (1400, 1600),
    "HA": (1600, 1800),
}


def esearch_ids(term: str, retmax: int) -> list[str]:
    handle = Entrez.esearch(db="nucleotide", term=term, retmax=retmax, usehistory="y")
    rec = Entrez.read(handle)
    handle.close()
    return list(rec["IdList"])


def efetch_gb(ids: list[str]) -> list:
    out = []
    for i in range(0, len(ids), 200):
        chunk = ids[i : i + 200]
        handle = Entrez.efetch(db="nucleotide", id=chunk, rettype="gb", retmode="text")
        out.extend(list(SeqIO.parse(handle, "genbank")))
        handle.close()
        time.sleep(0.35)
    return out


def extract_aa(rec) -> SeqRecord | None:
    """Translate first CDS, or whole seq if CDS missing."""
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        try:
            aa = feat.qualifiers.get("translation", [None])[0]
            if aa:
                acc = rec.id.split(".")[0]
                date = "2000-01-01"
                if rec.annotations.get("date"):
                    # GenBank date is like 01-JAN-2020 — keep year if possible
                    parts = str(rec.annotations["date"]).split("-")
                    if len(parts) == 3 and parts[2].isdigit():
                        date = f"{parts[2]}-01-01"
                return SeqRecord(Seq(aa), id=f"{acc},{date}", description="")
            nt = feat.extract(rec.seq)
            if len(nt) % 3:
                nt = nt[: len(nt) - (len(nt) % 3)]
            return SeqRecord(nt.translate(to_stop=False), id=rec.id, description="")
        except Exception:  # noqa: BLE001
            continue
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--subtype", required=True, choices=["H3N2", "H1N1"])
    ap.add_argument("--segment", required=True, choices=list(SEG_LEN))
    ap.add_argument("--retmax", type=int, default=8000)
    ap.add_argument("--host", default="Homo sapiens")
    args = ap.parse_args()

    lo, hi = SEG_LEN[args.segment]
    gene_id = f"{args.subtype.lower()}_{args.segment.lower()}"
    out_dir = ROOT / "data" / gene_id / "raw"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fa = out_dir / f"{gene_id}.fasta"

    term = (
        f'"Influenza A virus"[Organism] AND {args.subtype}[All Fields] '
        f"AND {args.segment}[Gene Name] AND \"{args.host}\"[Host] "
        f"AND {lo}:{hi}[SLEN] AND 1924:2026[PDAT]"
    )
    print(f"query: {term}")
    ids = esearch_ids(term, args.retmax)
    print(f"ids={len(ids)}")
    recs = efetch_gb(ids)
    written = 0
    seen: set[str] = set()
    with open(out_fa, "w") as fh:
        for rec in recs:
            aa = extract_aa(rec)
            if aa is None:
                continue
            rid = aa.id.split(",")[0]
            if rid in seen:
                continue
            seen.add(rid)
            SeqIO.write(aa, fh, "fasta")
            written += 1
    print(f"Wrote {out_fa} n={written}")


if __name__ == "__main__":
    main()
