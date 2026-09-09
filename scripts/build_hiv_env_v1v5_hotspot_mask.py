#!/usr/bin/env python3
"""
Build mut_hotspot_mask for HIV-1 Env HXB2 V1–V5 variable loops.

HXB2 Env AA (1-based, Los Alamos):
  V1 131–157, V2 158–196, V3 296–331, V4 385–418, V5 461–471

Mask is a BoolTensor[L] in HXB2 column frame (col = hxb2_pos - 1), same
convention as COVID Wuhan-frame / flu ORF-frame lit masks. Extra mut weight
at train time via --mut-hotspot-mask / --mut-hotspot-weight.

Output:
  {"mut_hotspot_mask": BoolTensor[L], "meta": {...}}
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parent.parent

# 1-based inclusive HXB2 Env coordinates
HXB2_V_REGIONS = {
    "V1": (131, 157),
    "V2": (158, 196),
    "V3": (296, 331),
    "V4": (385, 418),
    "V5": (461, 471),
}


def build_mask(max_seq_len: int) -> tuple[torch.Tensor, dict]:
    mask = torch.zeros(max_seq_len, dtype=torch.bool)
    covered = {}
    for name, (lo, hi) in HXB2_V_REGIONS.items():
        # 1-based inclusive → 0-based half-open
        a = max(0, lo - 1)
        b = min(max_seq_len, hi)
        if a < b:
            mask[a:b] = True
        covered[name] = {"hxb2_1based": [lo, hi], "cols_0based": [a, b]}
    meta = {
        "reference": "HXB2 Env (K03455)",
        "numbering": "HXB2 Env AA 1-based; col_0based = hxb2_pos - 1",
        "regions": covered,
        "max_seq_len": max_seq_len,
        "n_hot_sites": int(mask.sum().item()),
        "note": (
            "Assumes TreeSBM columns track ungapped HXB2-like Env frame "
            "(near-complete Env CDS, length ~856). Indel-rich V-loops may "
            "shift columns on MSA; mask is a lit region prior like COVID RBD band."
        ),
    }
    return mask, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--max-seq-len",
        type=int,
        default=900,
        help="Mask length / train max_seq_len (HIV Env ~856 AA)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results" / "hiv_env_mask" / "mut_hotspot_mask_hiv_env_v1v5.pt",
    )
    args = ap.parse_args()

    mask, meta = build_mask(args.max_seq_len)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    blob = {"mut_hotspot_mask": mask.cpu(), "meta": meta}
    torch.save(blob, args.out)
    meta_path = args.out.with_suffix(".meta.json")
    meta_path.write_text(json.dumps(meta, indent=2))
    print(f"Wrote {args.out}  hot={meta['n_hot_sites']}/{args.max_seq_len}")
    print(f"Wrote {meta_path}")


if __name__ == "__main__":
    main()
