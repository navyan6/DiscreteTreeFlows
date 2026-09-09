#!/usr/bin/env python3
"""
Feasibility probe: can we actually build time x population trees for a virus?

Splitting by "Ebola, 2016-2017, South Africa" only works if GenBank records
carry /collection_date and /geo_loc_name. Coverage varies enormously by virus,
so this samples records and reports the yield that survives the split keys.

The number that matters is not how many sequences a virus has, but how many
land in a single (country, time window) bucket with a usable date.
"""

from __future__ import annotations

import argparse
import collections
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request

from Bio import SeqIO

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# A deliberately mixed panel: different families, genome types, and very
# different surveillance intensities.
DEFAULT_PANEL = [
    (64320, "Zika virus", 9000, 12000),
    (12637, "Dengue virus", 9000, 12000),
    (11234, "Measles morbillivirus", 14000, 17000),
    (186538, "Zaire ebolavirus", 17000, 20000),
    (11676, "HIV-1", 8000, 10000),
    (11320, "Influenza A virus", 1500, 1900),
    (2697049, "SARS-CoV-2", 29000, 30500),
    (138948, "Enterovirus A", 6000, 8000),
    (11103, "Hepatitis C virus", 9000, 10000),
    (10298, "Human alphaherpesvirus 1", 140000, 160000),
    (114727, "Influenza A H1N1", 1500, 1900),
    (12092, "Hepatitis A virus", 7000, 8000),
]

YEAR = re.compile(r"(\d{4})")


def eget(endpoint: str, params: dict, timeout: int = 300) -> str:
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception:  # noqa: BLE001 - transient; back off and retry
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)
    return ""


def probe(tax_id: int, label: str, lo: int, hi: int, sample: int,
          sleep: float) -> dict:
    term = f"txid{tax_id}[Organism:exp] AND {lo}:{hi}[SLEN]"
    res = json.loads(eget("esearch.fcgi", {
        "db": "nuccore", "term": term, "retmax": sample, "retmode": "json",
    }))["esearchresult"]
    total = int(res["count"])
    ids = res["idlist"]
    time.sleep(sleep)

    if not ids:
        return {"tax_id": tax_id, "label": label, "total_in_range": total,
                "sampled": 0}

    gb = eget("efetch.fcgi", {"db": "nuccore", "id": ",".join(ids),
                              "rettype": "gb", "retmode": "text"})
    time.sleep(sleep)

    n = dated = placed = both = annotated = 0
    years: list[int] = []
    countries: collections.Counter = collections.Counter()
    buckets: collections.Counter = collections.Counter()

    for rec in SeqIO.parse(io.StringIO(gb), "genbank"):
        n += 1
        src = next((f for f in rec.features if f.type == "source"), None)
        q = src.qualifiers if src else {}
        date = (q.get("collection_date") or [None])[0]
        geo = (q.get("geo_loc_name") or q.get("country") or [None])[0]
        if any(f.type == "CDS" for f in rec.features):
            annotated += 1
        yr = None
        if date:
            m = YEAR.search(date)
            if m:
                yr = int(m.group(1))
                if 1900 < yr <= 2026:
                    dated += 1
                    years.append(yr)
                else:
                    yr = None
        if geo:
            placed += 1
            countries[geo.split(":")[0].strip()] += 1
        if yr and geo:
            both += 1
            # 2-year windows, the granularity of "2016-2017 in South Africa"
            buckets[(geo.split(":")[0].strip(), yr - yr % 2)] += 1

    frac = (lambda x: round(x / n, 3)) if n else (lambda x: 0.0)
    top = buckets.most_common(3)
    return {
        "tax_id": tax_id, "label": label, "total_in_range": total,
        "sampled": n,
        "frac_dated": frac(dated), "frac_placed": frac(placed),
        "frac_usable": frac(both), "frac_cds_annotated": frac(annotated),
        "year_range": [min(years), max(years)] if years else None,
        "frac_post_2024": round(sum(y > 2024 for y in years) / max(len(years), 1), 3),
        "n_countries": len(countries),
        "top_countries": countries.most_common(4),
        "top_buckets": [{"country": c, "window": f"{w}-{w+1}",
                         "frac_of_sample": round(k / n, 3),
                         "projected_seqs": int(round(total * both / n * (k / max(both, 1))))}
                        for (c, w), k in top],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", type=int, default=200)
    ap.add_argument("--sleep", type=float, default=0.5)
    ap.add_argument("--out", default="data/panviral/metadata_yield.json")
    args = ap.parse_args()

    rows = []
    hdr = (f"{'virus':<26}{'in-range':>10}{'dated':>7}{'geo':>7}"
           f"{'usable':>8}{'CDS':>6}  {'years':<12}{'largest (country, window)'}")
    print(hdr)
    print("-" * len(hdr))
    for tax_id, label, lo, hi in DEFAULT_PANEL:
        try:
            r = probe(tax_id, label, lo, hi, args.sample, args.sleep)
        except Exception as e:  # noqa: BLE001 - report and continue
            print(f"{label:<26} ERROR {e}")
            continue
        rows.append(r)
        if not r.get("sampled"):
            print(f"{label:<26}{r['total_in_range']:>10,}  (no records in range)")
            continue
        yr = f"{r['year_range'][0]}-{r['year_range'][1]}" if r["year_range"] else "-"
        b = r["top_buckets"][0] if r["top_buckets"] else None
        btxt = (f"{b['country'][:18]} {b['window']} ~{b['projected_seqs']:,}"
                if b else "-")
        print(f"{r['label']:<26}{r['total_in_range']:>10,}"
              f"{r['frac_dated']:>7.0%}{r['frac_placed']:>7.0%}"
              f"{r['frac_usable']:>8.0%}{r['frac_cds_annotated']:>6.0%}  "
              f"{yr:<12}{btxt}")
        sys.stdout.flush()

    with open(args.out, "w") as fh:
        json.dump(rows, fh, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
