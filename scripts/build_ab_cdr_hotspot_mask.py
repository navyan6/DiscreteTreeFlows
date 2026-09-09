#!/usr/bin/env python3
"""
Build a TreeSBM mut-hotspot mask that upweights IMGT CDR1/2/3 for antibodies.

Two modes:
  1) empirical (preferred): locate ANARCI CDR1/2/3 substrings in ungapped AA
     sequences and aggregate per-site CDR occupancy → bool mask [max_seq_len].
  2) stub: scale classic IMGT VH CDR ranges (27–38, 56–65, 105–117) on a
     128-ref frame to max_seq_len (same idea as scripts/ab_t6_metrics.py).

Output .pt dict matches pathogen lit masks:
  {"mut_hotspot_mask": BoolTensor[L], "meta": {...}}

Usage:
  python scripts/build_ab_cdr_hotspot_mask.py \\
    --fasta data/ab_t5_1m/sequences.fasta \\
    --anarci data/ab_t5_1m/anarci.tsv \\
    --max-seq-len 160 \\
    --out results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr.pt

  # stub only (no ANARCI):
  python scripts/build_ab_cdr_hotspot_mask.py --stub --max-seq-len 160 \\
    --out results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_stub.pt
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parent.parent

# IMGT CDR inclusive position labels (ANARCI numbering)
_IMGT_CDR_RANGES = [(27, 38), (56, 65), (105, 117)]


def iter_fasta(path: Path):
    header = None
    parts: list[str] = []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if header is not None:
            yield header, "".join(parts)


def find_cdr_spans(seq: str, cdr1: str, cdr2: str, cdr3: str) -> list[tuple[int, int]]:
    """Return [start, end) spans in ungapped seq for non-empty CDR strings."""
    spans: list[tuple[int, int]] = []
    cursor = 0
    for cdr in (cdr1, cdr2, cdr3):
        c = (cdr or "").upper().replace("-", "").replace(".", "")
        if not c:
            continue
        i = seq.find(c, cursor)
        if i < 0:
            i = seq.find(c)  # allow non-monotonic if numbering quirk
        if i < 0:
            continue
        spans.append((i, i + len(c)))
        cursor = i + len(c)
    return spans


def stub_mask(max_seq_len: int) -> torch.Tensor:
    scale = max_seq_len / 128.0
    mask = torch.zeros(max_seq_len, dtype=torch.bool)
    for a, b in _IMGT_CDR_RANGES:
        lo = max(0, int(a * scale) - 1)
        hi = min(max_seq_len, int(b * scale))
        mask[lo:hi] = True
    return mask


def empirical_mask(
    fasta: Path,
    anarci: Path,
    max_seq_len: int,
    min_occupancy: float,
    max_seqs: int | None,
) -> tuple[torch.Tensor, dict]:
    by_id: dict[str, dict] = {}
    with anarci.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sid = row.get("sequence_id") or ""
            if sid:
                by_id[sid] = row

    occ = torch.zeros(max_seq_len, dtype=torch.float64)
    n_used = 0
    n_miss = 0
    n_no_cdr = 0
    cdr_len_hist: Counter[str] = Counter()

    for i, (sid, seq) in enumerate(iter_fasta(fasta)):
        if max_seqs is not None and i >= max_seqs:
            break
        row = by_id.get(sid)
        if row is None:
            n_miss += 1
            continue
        seq_u = seq.upper().replace("-", "").replace(".", "")
        spans = find_cdr_spans(
            seq_u,
            row.get("cdr1", ""),
            row.get("cdr2", ""),
            row.get("cdr3", ""),
        )
        if not spans:
            n_no_cdr += 1
            continue
        for a, b in spans:
            a2 = max(0, min(max_seq_len, a))
            b2 = max(0, min(max_seq_len, b))
            if b2 > a2:
                occ[a2:b2] += 1.0
                cdr_len_hist[str(b - a)] += 1
        n_used += 1

    if n_used == 0:
        raise SystemExit("No sequences with locatable CDRs — cannot build empirical mask")

    frac = occ / float(n_used)
    mask = frac >= float(min_occupancy)
    meta = {
        "mode": "empirical_anarci_substring",
        "n_used": n_used,
        "n_miss_anarci": n_miss,
        "n_no_cdr_span": n_no_cdr,
        "min_occupancy": min_occupancy,
        "n_hot_sites": int(mask.sum().item()),
        "frac_hot": float(mask.float().mean().item()),
        "mean_occupancy_hot": float(frac[mask].mean().item()) if mask.any() else 0.0,
        "cdr_span_len_top": dict(cdr_len_hist.most_common(15)),
        "imgt_ranges_ref": _IMGT_CDR_RANGES,
    }
    return mask, meta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=Path, default=None)
    ap.add_argument("--anarci", type=Path, default=None)
    ap.add_argument("--max-seq-len", type=int, default=160)
    ap.add_argument("--min-occupancy", type=float, default=0.25,
                    help="Empirical: site is CDR-hot if ≥ this fraction of seqs mark it CDR")
    ap.add_argument("--max-seqs", type=int, default=None, help="Cap sequences scanned (empirical)")
    ap.add_argument("--stub", action="store_true", help="Use scaled IMGT stub (no ANARCI)")
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results" / "ab_cdr_mask" / "mut_hotspot_mask_imgt_cdr.pt",
    )
    args = ap.parse_args()

    if args.stub:
        mask = stub_mask(args.max_seq_len)
        meta = {
            "mode": "imgt_stub_scaled",
            "max_seq_len": args.max_seq_len,
            "imgt_ranges_ref128": _IMGT_CDR_RANGES,
            "n_hot_sites": int(mask.sum().item()),
            "frac_hot": float(mask.float().mean().item()),
        }
    else:
        if args.fasta is None or args.anarci is None:
            raise SystemExit("Empirical mode requires --fasta and --anarci (or pass --stub)")
        mask, meta = empirical_mask(
            args.fasta, args.anarci, args.max_seq_len, args.min_occupancy, args.max_seqs
        )
        meta["max_seq_len"] = args.max_seq_len
        meta["fasta"] = str(args.fasta)
        meta["anarci"] = str(args.anarci)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    blob = {"mut_hotspot_mask": mask.cpu(), "meta": meta}
    torch.save(blob, args.out)
    meta_path = args.out.with_suffix(".meta.json")
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    print(f"Wrote {args.out}  hot_sites={int(mask.sum().item())}/{args.max_seq_len}")


if __name__ == "__main__":
    main()
