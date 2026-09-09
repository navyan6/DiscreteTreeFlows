#!/usr/bin/env python3
"""
Ingest Nextstrain filovirus metadata via the Charon auspice API.

Public Nextstrain policy: metadata + accessions are downloadable; whole-genome
sequences are fetched from NCBI by accession (see bdbv_extract_l.py).

Merges outbreak labels into data/bdbv/manifest.json (NCBI rows preferred on conflict).

Usage:
  python scripts/download_ebolavirus_nextstrain.py
  python scripts/download_ebolavirus_nextstrain.py --build ebola/bdbv
  python scripts/download_ebolavirus_nextstrain.py --dry-run
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.bdbv_common import (
    DATA_ROOT,
    dedupe_manifest,
    enrich_manifest_outbreaks,
    infer_outbreak_id,
    load_manifest,
    normalize_country,
    parse_date,
    save_manifest,
    slug_outbreak_id,
    write_data_audit,
)

CHARON_URL = "https://nextstrain.org/charon/getDataset?prefix={prefix}"

# (charon prefix, default species or "auto", optional fixed outbreak_id slug)
NEXTSTRAIN_BUILDS: list[tuple[str, str, str | None]] = [
    ("ebola/ebov-2013", "ebov", "ebov_wa_2013_2016"),
    ("community/inrb-drc/ebola-nord-kivu", "ebov", "ebov_nordkivu_2018_2020"),
    ("ebola/all-outbreaks", "auto", None),
    ("ebola/bdbv", "bdbv", None),
    ("ebola/bdbv-drc-uganda-2026", "bdbv", "bdbv_2026"),
    ("ebola/sudv", "sudv", None),
]

ACCESSION_RE = re.compile(r"^[A-Z]{1,2}\d+[A-Z]?\d*$")


def _resolve_accession(attrs: dict[str, Any]) -> str:
    """Prefer INSDC/GenBank accessions (NCBI efetch); ignore Pathoplexus-only IDs."""
    for key in ("INSDC_accession", "genbank_accession", "accession"):
        raw = str(_attr_value(attrs, key) or "").strip()
        if not raw:
            continue
        acc = raw.split(".")[0].split("|")[0]
        if ACCESSION_RE.match(acc):
            return acc
    strain = str(_attr_value(attrs, "strain") or "")
    for token in re.split(r"[/|,\s]+", strain):
        tok = token.split(".")[0]
        if ACCESSION_RE.match(tok):
            return tok
    return ""


def _resolve_date(attrs: dict[str, Any]) -> tuple[str, int | None]:
    date_raw = str(_attr_value(attrs, "date") or "").strip()
    if date_raw:
        augur, _, year = parse_date(date_raw)
        return augur, year
    return _num_date_to_augur(_attr_value(attrs, "num_date"))


def _attr_value(node_attrs: dict[str, Any], key: str) -> Any:
    val = node_attrs.get(key)
    if isinstance(val, dict):
        return val.get("value", "")
    return val


def _num_date_to_augur(num_date: Any) -> tuple[str, int | None]:
    if num_date is None or num_date == "":
        return "2000-XX-XX", None
    if isinstance(num_date, dict):
        num_date = num_date.get("value", "")
    try:
        yf = float(num_date)
        year = int(yf)
        return f"{year:04d}-XX-XX", year
    except (TypeError, ValueError):
        augur, _, year = parse_date(str(num_date))
        return augur, year


def _species_from_text(*parts: str) -> str:
    txt = " ".join(p for p in parts if p).lower()
    if "bundibugyo" in txt or "bdbv" in txt:
        return "bdbv"
    if "reston" in txt:
        return "restv"
    if "marburg" in txt or "ravn" in txt:
        return "marv"
    if "sudan" in txt and "ebola" in txt:
        return "sudv"
    if "tai" in txt and "forest" in txt:
        return "tafv"
    return "ebov"


def _outbreak_from_build(
    build_prefix: str,
    fixed_outbreak: str | None,
    species: str,
    division: str,
    country: str,
    year: int | None,
) -> str:
    if fixed_outbreak:
        return fixed_outbreak
    if division and division.lower() not in ("?", "unknown", "unassigned"):
        return slug_outbreak_id(f"{species}_{division}")
    tail = build_prefix.split("/")[-1].replace("-", "_")
    if tail and tail not in ("ebola", "all_outbreaks"):
        return slug_outbreak_id(f"{species}_{tail}")
    return infer_outbreak_id(
        {"species": species, "country": country, "year": year, "outbreak": division}
    )


def fetch_charon_dataset(prefix: str, timeout: int = 120) -> dict[str, Any]:
    url = CHARON_URL.format(prefix=prefix)
    req = urllib.request.Request(url, headers={"User-Agent": "DiscreteTreeFlows/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def iter_leaves(node: dict[str, Any]) -> list[dict[str, Any]]:
    children = node.get("children") or []
    if not children:
        return [node]
    out: list[dict[str, Any]] = []
    for child in children:
        out.extend(iter_leaves(child))
    return out


def ingest_build(
    prefix: str,
    default_species: str,
    fixed_outbreak: str | None,
) -> list[dict[str, Any]]:
    print(f"Charon {prefix} ...", flush=True)
    try:
        payload = fetch_charon_dataset(prefix)
    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError) as e:
        print(f"  FAIL fetch: {e}")
        return []

    tree = payload.get("tree")
    if not tree:
        print("  FAIL: no tree in payload")
        return []

    leaves = iter_leaves(tree)
    records: list[dict[str, Any]] = []
    skipped = 0
    for leaf in leaves:
        attrs = leaf.get("node_attrs") or {}
        acc = _resolve_accession(attrs)
        if not acc:
            skipped += 1
            continue

        country_raw = str(_attr_value(attrs, "country") or "")
        division = str(_attr_value(attrs, "division") or "")
        country = normalize_country(country_raw or division or "unknown")
        augur_date, year = _resolve_date(attrs)

        species = default_species
        if species == "auto":
            species = _species_from_text(prefix, division, country, leaf.get("name", ""))

        outbreak_label = str(_attr_value(attrs, "outbreak") or division or "")
        outbreak_id = _outbreak_from_build(
            prefix, fixed_outbreak, species, outbreak_label or division, country, year
        )
        strain = str(_attr_value(attrs, "strain") or leaf.get("name") or acc)
        ppx = str(_attr_value(attrs, "PPX_accession") or "")

        records.append(
            {
                "accession": acc,
                "species": species,
                "date": augur_date,
                "year": year,
                "country": country,
                "source": "nextstrain",
                "host": "",
                "strain": strain,
                "division": division,
                "outbreak": outbreak_label or division or prefix.split("/")[-1],
                "outbreak_id": outbreak_id,
                "nextstrain_build": prefix,
                "pathoplexus_id": ppx,
            }
        )

    print(f"  leaves={len(leaves)}  with_accession={len(records)}  skipped={skipped}")
    return records


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--build",
        action="append",
        default=[],
        help="Charon prefix (repeatable). Default: all NEXTSTRAIN_BUILDS.",
    )
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    builds = NEXTSTRAIN_BUILDS
    if args.build:
        wanted = set(args.build)
        builds = [b for b in NEXTSTRAIN_BUILDS if b[0] in wanted]
        missing = wanted - {b[0] for b in builds}
        for prefix in sorted(missing):
            builds.append((prefix, "auto", None))

    new_records: list[dict] = []
    for prefix, default_sp, fixed_ob in builds:
        new_records.extend(ingest_build(prefix, default_sp, fixed_ob))

    # Dedupe within nextstrain batch (same acc may appear in multiple builds)
    by_acc: dict[str, dict] = {}
    for r in new_records:
        acc = r["accession"]
        prev = by_acc.get(acc)
        if prev is None:
            by_acc[acc] = dict(r)
            continue
        for key in ("outbreak", "outbreak_id", "division", "strain", "nextstrain_build"):
            if not prev.get(key) and r.get(key):
                prev[key] = r[key]
    new_records = list(by_acc.values())

    existing = load_manifest()
    before = len(existing)
    merged = dedupe_manifest(existing + new_records)
    merged = enrich_manifest_outbreaks(merged)
    added = len(merged) - before

    print(
        f"Nextstrain metadata: {len(new_records)} unique accessions "
        f"({len(new_records) - added} overlapped NCBI); manifest {before} -> {len(merged)} (+{added})"
    )

    if args.dry_run:
        print("Dry run — manifest not written.")
        return

    save_manifest(merged)
    write_data_audit(merged)
    audit_path = DATA_ROOT / "nextstrain_charon_audit.json"
    audit_path.write_text(
        json.dumps(
            {
                "builds": [b[0] for b in builds],
                "nextstrain_accessions": len(new_records),
                "manifest_total": len(merged),
                "new_accessions": added,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"Wrote {DATA_ROOT / 'manifest.json'} and {audit_path}")


if __name__ == "__main__":
    main()
