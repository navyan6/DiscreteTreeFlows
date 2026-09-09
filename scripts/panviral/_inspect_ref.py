#!/usr/bin/env python3
"""Scratch helper: dump the annotated features of a RefSeq accession."""
import collections
import io
import sys
import urllib.parse
import urllib.request

from Bio import SeqIO

acc = sys.argv[1]
url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(
    {"db": "nuccore", "id": acc, "rettype": "gb", "retmode": "text"}
)
gb = urllib.request.urlopen(url, timeout=120).read().decode()
rec = next(SeqIO.parse(io.StringIO(gb), "genbank"))
print(acc, len(rec.seq), "nt")
print("feature types:", dict(collections.Counter(f.type for f in rec.features)))
for f in rec.features:
    if f.type in ("mat_peptide", "CDS", "sig_peptide"):
        prod = f.qualifiers.get("product", ["?"])[0]
        print(f"  {f.type:12s} {int(f.location.start):>6}-{int(f.location.end):<6} {prod[:60]}")
