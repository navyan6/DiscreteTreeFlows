"""Shared helpers for Bundibugyo / filovirus L-protein TreeSBM pipeline."""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = ROOT / "data" / "bdbv"

# Reference L CDS coordinates (0-based half-open on NC_014373.1 / BDBV)
REF_BDBV_ACCESSION = "NC_014373.1"
REF_BDBV_L_NT = (11565, 18198)  # 2211 codons
REF_FJ217161 = "FJ217161.1"

# Track A gold holdout — must never appear in train/val (eval only).
GOLDEN_TEST_OUTBREAK = "bdbv_2026"
GOLDEN_TEST_OUTBREAKS = frozenset({GOLDEN_TEST_OUTBREAK})

# EBOV Makona KM034562 L CDS (typical; verified via GenBank feature)
REF_EBOV_ACCESSION = "KM034562.1"
REF_EBOV_L_NT = (11508, 18140)  # ~2211 aa

SPECIES: dict[str, dict[str, Any]] = {
    "ebov": {
        "label": "EBOV",
        "organism": "Orthoebolavirus zairense",
        "taxid": 186538,
        "ref_accession": REF_EBOV_ACCESSION,
        "l_nt": (11579, 18219),
    },
    "sudv": {
        "label": "SUDV",
        "organism": "Orthoebolavirus sudanense",
        "taxid": 186540,
        "ref_accession": "NC_006432.1",
        "l_nt": (11508, 18140),
    },
    "bdbv": {
        "label": "BDBV",
        "organism": "Orthoebolavirus bundibugyoense",
        "taxid": 565995,
        "ref_accession": REF_BDBV_ACCESSION,
        "l_nt": (11565, 18198),
    },
    "tafv": {
        "label": "TAFV",
        "organism": "Orthoebolavirus taiense",
        "taxid": 186539,
        "ref_accession": "NC_014372.1",
        "l_nt": (11508, 18140),
    },
    "restv": {
        "label": "RESTV",
        "organism": "Orthoebolavirus restonense",
        "taxid": 186541,
        "ref_accession": "NC_004161.1",
        "l_nt": (11508, 18140),
    },
    "marv": {
        "label": "MARV",
        "organism": "Marburgvirus marburgvirus",
        "taxid": 11269,
        "ref_accession": "NC_001608.1",
        "l_nt": (11552, 18185),
    },
}

WEST_AFRICA_COUNTRIES = {
    "guinea",
    "sierra leone",
    "liberia",
    "senegal",
    "mali",
    "nigeria",
}

EAST_AFRICA_COUNTRIES = {
    "uganda",
    "sudan",
    "south sudan",
    "democratic republic of the congo",
    "drc",
    "congo",
}


@dataclass
class GenomeRecord:
    accession: str
    species: str
    date: str
    country: str
    source: str
    length: int
    host: str = ""

    def header_suffix(self) -> str:
        return f"{self.date},{self.country},{self.species},{self.source}"


def parse_date(raw: str) -> tuple[str, tuple[int, int, int], int | None]:
    """Parse collection date -> (augur str, sort key, year)."""
    raw = (raw or "").strip()
    if not raw:
        return "2000-XX-XX", (2000, 7, 15), None
    parts = re.split(r"[-/]", raw)
    try:
        y = int(parts[0])
    except ValueError:
        return "2000-XX-XX", (2000, 7, 15), None
    if len(parts) == 1:
        return f"{y:04d}-XX-XX", (y, 7, 15), y
    if len(parts) == 2:
        if parts[1].upper() == "XX":
            return f"{y:04d}-XX-XX", (y, 7, 15), y
        m = int(parts[1])
        return f"{y:04d}-{m:02d}-XX", (y, m, 15), y
    m_part, d_part = parts[1], parts[2]
    if m_part.upper() == "XX":
        return f"{y:04d}-XX-XX", (y, 7, 15), y
    m = int(m_part)
    if d_part.upper() == "XX":
        return f"{y:04d}-{m:02d}-XX", (y, m, 15), y
    d = int(d_part)
    return f"{y:04d}-{m:02d}-{d:02d}", (y, m, d), y


def normalize_country(raw: str) -> str:
    s = (raw or "").strip().lower()
    s = re.sub(r"\s+", " ", s)
    if s in ("drc", "dem. rep. congo", "democratic republic of the congo"):
        return "democratic republic of the congo"
    return s


def is_east_africa(country: str) -> bool:
    c = normalize_country(country)
    return any(x in c for x in EAST_AFRICA_COUNTRIES)


def manifest_path() -> Path:
    return DATA_ROOT / "manifest.json"


def load_manifest() -> list[dict[str, Any]]:
    p = manifest_path()
    if not p.is_file():
        return []
    return json.loads(p.read_text())


def save_manifest(records: list[dict[str, Any]]) -> None:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    manifest_path().write_text(json.dumps(records, indent=2) + "\n")


def window_config_path() -> Path:
    return ROOT / "results" / "bdbv_l_conservation" / "window_config.json"


def load_window_config() -> dict[str, Any]:
    p = window_config_path()
    if not p.is_file():
        return {"max_seq_len": 900, "aa_start": 0, "aa_end": 900}
    return json.loads(p.read_text())


def slug_outbreak_id(raw: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9]+", "_", (raw or "").strip()).strip("_").lower()
    return s[:96] or "unknown"


def infer_outbreak_id(record: dict[str, Any]) -> str:
    """Assign one outbreak key per genome (Nextstrain label or heuristic)."""
    for key in ("outbreak_id", "outbreak"):
        val = (record.get(key) or "").strip()
        if val and val.lower() not in ("?", "unknown", "na", "n/a"):
            sp = record.get("species", "unknown")
            return slug_outbreak_id(f"{sp}_{val}")

    sp = record.get("species", "unknown")
    y = record.get("year")
    c = normalize_country(record.get("country", "unknown"))

    if sp == "bdbv":
        if y is not None and y >= 2026:
            return "bdbv_2026"
        if y is not None and y <= 2008:
            return "bdbv_2007"
        if y is not None and 2011 <= y <= 2013:
            return "bdbv_2012"

    if sp == "ebov" and y is not None:
        if 2013 <= y <= 2016 and (
            any(w in c for w in WEST_AFRICA_COUNTRIES) or c == "unknown"
        ):
            return "ebov_wa_2013_2016"
        if 2018 <= y <= 2020 and "congo" in c:
            return "ebov_nordkivu_2018_2020"
        if y == 1976:
            return "ebov_1976_yambuku"
        if 2014 <= y <= 2015 and "congo" in c:
            return "ebov_2014_drc"

    if sp == "sudv" and y is not None and 2010 <= y <= 2012:
        return "sudv_2011_uganda"

    if sp == "marv" and y is not None:
        if 2004 <= y <= 2005 and ("angola" in c or c == "unknown"):
            return "marv_angola_2005"
        if 1967 <= y <= 1968:
            return "marv_1967_marburg"

    if y is not None:
        return slug_outbreak_id(f"{sp}_{c}_{y}")
    return slug_outbreak_id(f"{sp}_{c}_unknown")


def enrich_manifest_outbreaks(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Fill outbreak_id on every record; merge labels from all sources on dedupe."""
    out: list[dict[str, Any]] = []
    for r in records:
        rec = dict(r)
        if not rec.get("outbreak_id"):
            rec["outbreak_id"] = infer_outbreak_id(rec)
        out.append(rec)
    return out


def write_data_audit(records: list[dict[str, Any]], out: Path | None = None) -> None:
    out = out or DATA_ROOT / "DATA_AUDIT.md"
    records = enrich_manifest_outbreaks(records)
    by_sp: dict[str, int] = {}
    by_year: dict[str, int] = {}
    by_src: dict[str, int] = {}
    by_outbreak: dict[str, int] = {}
    for r in records:
        sp = r.get("species", "?")
        by_sp[sp] = by_sp.get(sp, 0) + 1
        by_src[r.get("source", "?")] = by_src.get(r.get("source", "?"), 0) + 1
        oid = r.get("outbreak_id") or infer_outbreak_id(r)
        by_outbreak[oid] = by_outbreak.get(oid, 0) + 1
        y = r.get("year")
        if y is not None:
            by_year[str(y)] = by_year.get(str(y), 0) + 1
    lines = [
        "# BDBV / filovirus data audit",
        "",
        f"Total genomes in manifest: **{len(records)}**",
        "",
        "## By species",
        "",
        "| species | N |",
        "|---|---:|",
    ]
    for sp in sorted(by_sp):
        lines.append(f"| {sp} | {by_sp[sp]} |")
    lines += ["", "## By source", "", "| source | N |", "|---|---:|"]
    for src in sorted(by_src):
        lines.append(f"| {src} | {by_src[src]} |")
    lines += ["", "## By year", "", "| year | N |", "|---|---:|"]
    for y in sorted(by_year, key=lambda x: int(x) if x.isdigit() else 0):
        lines.append(f"| {y} | {by_year[y]} |")
    lines += ["", "## By outbreak (top 40)", "", "| outbreak_id | N |", "|---|---:|"]
    for oid, n in sorted(by_outbreak.items(), key=lambda x: -x[1])[:40]:
        lines.append(f"| {oid} | {n} |")
    lines.append("")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")


def extract_l_nt_from_genbank(rec) -> str | None:
    """Extract L-gene NT CDS from a GenBank record (gene=L feature)."""
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        gene = feat.qualifiers.get("gene", [""])[0]
        if gene != "L":
            continue
        region = str(feat.extract(rec.seq)).replace("-", "").upper()
        if len(region) < 5500 or len(region) % 3 != 0:
            continue
        from Bio.Seq import Seq

        tr = str(Seq(region).translate())
        if "*" in tr[:-1]:
            continue
        return region
    return None


def parse_fasta_header(rec) -> tuple[str, str, str, str, str]:
    """Parse >ACC date,country,species,source or legacy >ACC,date,... comma-only."""
    acc = rec.id.split()[0].strip()
    if rec.description and rec.description.strip():
        desc = rec.description.strip()
        if desc.startswith(acc):
            desc = desc[len(acc) :].strip()
        parts = [acc] + [p.strip() for p in desc.split(",")]
    else:
        parts = [p.strip() for p in rec.id.split(",")]
    date = parts[1] if len(parts) > 1 else "2000-XX-XX"
    country = parts[2] if len(parts) > 2 else "unknown"
    species = parts[3] if len(parts) > 3 else "bdbv"
    source = parts[4] if len(parts) > 4 else "ncbi"
    return acc, date, country, species, source


def format_fasta_header(acc: str, date: str, country: str, species: str, source: str) -> str:
    return f">{acc} {date},{country},{species},{source}"


def dedupe_manifest(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Prefer ncbi > nextstrain > pathoplexus on duplicate accession; merge outbreak metadata."""
    rank = {"ncbi": 0, "nextstrain": 1, "pathoplexus": 2}
    best: dict[str, dict[str, Any]] = {}
    extras: dict[str, list[dict[str, Any]]] = {}
    for r in records:
        acc = r["accession"].split(".")[0]
        extras.setdefault(acc, []).append(r)
        prev = best.get(acc)
        if prev is None or rank.get(r.get("source", ""), 9) < rank.get(
            prev.get("source", ""), 9
        ):
            best[acc] = dict(r)

    for acc, rec in best.items():
        for other in extras.get(acc, []):
            for key in (
                "outbreak",
                "outbreak_id",
                "strain",
                "division",
                "pathoplexus_id",
                "golden_test",
            ):
                if not rec.get(key) and other.get(key):
                    rec[key] = other[key]
        # Pathoplexus / Nextstrain 2026 labels win over NCBI heuristics on dedupe.
        for other in extras.get(acc, []):
            oid = other.get("outbreak_id") or ""
            ob = (other.get("outbreak") or "").lower()
            if oid in GOLDEN_TEST_OUTBREAKS or ob == "bdbv-2026":
                rec["outbreak_id"] = GOLDEN_TEST_OUTBREAK
                rec["outbreak"] = other.get("outbreak") or "Bdbv-2026"
                rec["golden_test"] = True
                if other.get("pathoplexus_id"):
                    rec["pathoplexus_id"] = other["pathoplexus_id"]
                break
        if not rec.get("outbreak_id"):
            if rec.get("outbreak"):
                rec["outbreak_id"] = slug_outbreak_id(
                    f"{rec.get('species', 'unknown')}_{rec['outbreak']}"
                )
            else:
                rec["outbreak_id"] = infer_outbreak_id(rec)
    return sorted(enrich_manifest_outbreaks(list(best.values())), key=lambda x: x["accession"])
