#!/usr/bin/env python3
"""
Table 5 (Ab) — subsample ~500K AA sequences from data/Homo_sapiens.fasta.

Prefer PairedNGS / BioProject repertoire dumps (SRR = donor/person proxy),
heavy chain, then fill with other sources if needed. Writes a clean FASTA +
metadata TSV for ANARCI / clone clustering.

Does NOT run ANARCI, define clones, or build trees.

Usage:
  python scripts/prepare_ab_t5_subsample.py \\
    --input data/Homo_sapiens.fasta \\
    --out-dir data/ab_t5_500k \\
    --max-seqs 500000
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import importlib.util

_spec = importlib.util.spec_from_file_location(
    "prepare_ab_fasta", ROOT / "scripts" / "prepare_ab_fasta.py"
)
_pab = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_pab)
_BIOPROJECT_RE = _pab._BIOPROJECT_RE
_CLONE_TAG_RE = _pab._CLONE_TAG_RE
_SRR_RE = _pab._SRR_RE
_VGENE_RE = _pab._VGENE_RE
infer_chain = _pab.infer_chain
iter_fasta = _pab.iter_fasta
sanitize_aa = _pab.sanitize_aa

_PAIRED_RE = re.compile(r"PairedNGS", re.I)
_PATIENT_RE = re.compile(
    r"(Human Patient|Patient[_-]?\w+|Donor[_-]?\w+|Subject[_-]?\w+)", re.I
)
# PairedNGS block: PRJ…|SRR…|human|BARCODE|IGHV…|IGK…
_PAIRED_BLOCK_RE = re.compile(
    r"(PRJ[NE][A-Z]?\d+)\|(SRR\d+)\|([^|]*)\|([^|]*)\|([^|]*)",
    re.I,
)


def parse_donor(header: str) -> dict:
    """donor_id prefers SRR (per-person/run), else bioproject, else patient token."""
    h = header.strip()
    srr_m = _SRR_RE.search(h)
    bp_m = _BIOPROJECT_RE.search(h)
    pat_m = _PATIENT_RE.search(h)
    paired = bool(_PAIRED_RE.search(h) or (bp_m and srr_m))
    barcode = ""
    germline_h = ""
    block = _PAIRED_BLOCK_RE.search(h)
    if block:
        barcode = block.group(4).strip()
        germline_h = block.group(5).strip()

    if srr_m:
        donor_id = srr_m.group(1)
        donor_source = "srr"
    elif pat_m and pat_m.group(1).lower() not in {"human patient", "patient", "person"}:
        donor_id = pat_m.group(1)
        donor_source = "patient_token"
    elif bp_m:
        donor_id = bp_m.group(1)
        donor_source = "bioproject"
    else:
        donor_id = "unknown"
        donor_source = "unknown"

    v_gene = ""
    vm = _VGENE_RE.search(h)
    if vm:
        v_gene = vm.group(1).upper().replace("_", "-")
    if not v_gene and germline_h:
        gm = _VGENE_RE.search(germline_h.replace("_", "-"))
        if gm:
            v_gene = gm.group(1).upper().replace("_", "-")

    chain = infer_chain(h, v_gene)
    clone_m = _CLONE_TAG_RE.search(h)
    primary = h.split("|||", 1)[0].strip()
    seq_id = primary.split("|", 1)[0].split(";", 1)[0].strip() or "unknown"
    seq_id = re.sub(r"[^\w.\-]", "_", seq_id)[:120]

    return {
        "seq_id": seq_id,
        "donor_id": donor_id,
        "donor_source": donor_source,
        "bioproject": bp_m.group(1) if bp_m else "",
        "srr": srr_m.group(1) if srr_m else "",
        "barcode": barcode,
        "v_gene": v_gene,
        "chain": chain,
        "paired_ngs": int(paired),
        "clone_tag": clone_m.group(1) if clone_m else "",
        "header_short": primary[:200],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, default=ROOT / "data" / "Homo_sapiens.fasta")
    ap.add_argument("--out-dir", type=Path, default=ROOT / "data" / "ab_t5_500k")
    ap.add_argument("--max-seqs", type=int, default=500_000)
    ap.add_argument("--min-len", type=int, default=80)
    ap.add_argument("--max-len", type=int, default=160)
    ap.add_argument("--chains", default="heavy", help="Comma list; default heavy-only for T5")
    ap.add_argument("--prefer-paired", action="store_true", default=True)
    ap.add_argument("--no-prefer-paired", action="store_false", dest="prefer_paired")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--progress-every", type=int, default=200_000)
    args = ap.parse_args()

    allow = {c.strip().lower() for c in args.chains.split(",") if c.strip()}
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    # Reservoir for non-preferred fill; preferred stream kept greedily until full.
    preferred: list[tuple[dict, str]] = []
    other: list[tuple[dict, str]] = []
    n_read = n_filt = 0

    print(f"Reading {args.input} → reservoir up to {args.max_seqs:,}", flush=True)
    for header, raw in iter_fasta(args.input):
        n_read += 1
        if args.progress_every and n_read % args.progress_every == 0:
            print(
                f"  … read={n_read:,} pref={len(preferred):,} other={len(other):,}",
                flush=True,
            )
        aa = sanitize_aa(raw)
        L = len(aa)
        if L < args.min_len or L > args.max_len:
            n_filt += 1
            continue
        # reject obvious NT leftovers
        if set(aa) <= set("ACGTUN"):
            n_filt += 1
            continue
        meta = parse_donor(header)
        if meta["chain"] not in allow:
            n_filt += 1
            continue
        item = (meta, aa)
        is_pref = args.prefer_paired and meta["paired_ngs"] and meta["srr"]
        if is_pref:
            if len(preferred) < args.max_seqs:
                preferred.append(item)
            elif rng.random() < args.max_seqs / (len(preferred) + 1):
                # mild shuffle of preferred once full (keep diversity across file)
                preferred[rng.randrange(args.max_seqs)] = item
        else:
            # classic reservoir for fill
            if len(other) < args.max_seqs:
                other.append(item)
            else:
                j = rng.randrange(n_read)
                if j < args.max_seqs:
                    other[j % len(other)] = item

        if len(preferred) >= args.max_seqs:
            # still scan a bit more for diversity? stop early for <12h budget
            # Prefer stopping once preferred is full.
            print(f"Preferred pool full at read={n_read:,}; stopping early.", flush=True)
            break

    selected = preferred[: args.max_seqs]
    if len(selected) < args.max_seqs:
        need = args.max_seqs - len(selected)
        rng.shuffle(other)
        selected.extend(other[:need])

    fasta_path = out / "sequences.fasta"
    meta_path = out / "metadata.tsv"
    donor_counts: Counter[str] = Counter()
    chain_counts: Counter[str] = Counter()
    source_pref = 0

    with fasta_path.open("w") as ff, meta_path.open("w", newline="") as mf:
        w = csv.DictWriter(
            mf,
            fieldnames=[
                "uniq_id",
                "seq_id",
                "donor_id",
                "donor_source",
                "bioproject",
                "srr",
                "barcode",
                "v_gene",
                "chain",
                "paired_ngs",
                "clone_tag",
                "len",
                "header_short",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for i, (meta, aa) in enumerate(selected):
            uniq = f"t5ab{i:07d}_{meta['seq_id']}"[:80]
            ff.write(f">{uniq}\n")
            for k in range(0, len(aa), 80):
                ff.write(aa[k : k + 80] + "\n")
            row = {**meta, "uniq_id": uniq, "len": len(aa)}
            w.writerow(row)
            donor_counts[meta["donor_id"]] += 1
            chain_counts[meta["chain"]] += 1
            source_pref += int(meta["paired_ngs"])

    summary = {
        "input": str(args.input),
        "out_dir": str(out),
        "n_read": n_read,
        "n_filt_len_or_chain": n_filt,
        "n_kept": len(selected),
        "n_preferred_paired": source_pref,
        "n_donors": len(donor_counts),
        "chain_counts": dict(chain_counts),
        "donor_top": dict(donor_counts.most_common(20)),
        "max_seqs": args.max_seqs,
        "note": (
            "AA subsample for Table 5. Next: ANARCI annotate → clone define "
            "(donor+V+CDR3) → trees → clonal-lineage holdout."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in (
        "n_read", "n_kept", "n_preferred_paired", "n_donors", "chain_counts"
    )}, indent=2))
    print(f"Wrote {fasta_path} and {meta_path}", flush=True)


if __name__ == "__main__":
    main()
