#!/usr/bin/env python3
"""
Build literature-based L polymerase mut-hotspot mask for TreeSBM training.

Maps EBOV Makona reference positions to window indices via MSA, writes:
  results/bdbv_l_mask/mut_hotspot_mask_lit.pt
  results/bdbv_l_mask/mut_hotspot_mask_lit.meta.json
  results/bdbv_l_mask/mut_hotspot_mask_eval_variants.json

Literature regions (EBOV Makona 1-based aa):
  D759 ±5 adaptive hotspot
  RdRp domain ~540-850 (cluster peaks from MSA entropy)
  N-term 1-380 interaction hub
  CTD last 250 aa
  GDNQ 741-744 excluded (conservation core)
  BDBV 2026 K1738, Q1770
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import torch
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq

from scripts.bdbv_common import (
    DATA_ROOT,
    REF_EBOV_ACCESSION,
    REF_FJ217161,
    load_window_config,
    window_config_path,
)

GDNQ_EXCLUDE = (741, 744)


def _load_ebov_ref_aa() -> str:
    ref_path = DATA_ROOT / "references" / f"{REF_EBOV_ACCESSION}.fasta"
    if ref_path.is_file():
        return str(SeqIO.read(ref_path, "fasta").seq)
    # fallback: translate L from genome slice
    from scripts.bdbv_common import SPECIES

    lo, hi = SPECIES["ebov"]["l_nt"]
    rec = SeqIO.read(ref_path, "fasta") if ref_path.is_file() else None
    if rec is None:
        raise FileNotFoundError(f"Need reference {ref_path}")
    return str(Seq(str(rec.seq)[lo:hi]).translate())


def _ref_pos_to_msa_col(msa, ref_row_idx: int, aa_1based: int) -> int | None:
    seq = str(msa[ref_row_idx].seq)
    pos = 0
    for i, c in enumerate(seq):
        if c == "-":
            continue
        pos += 1
        if pos == aa_1based:
            return i
    return None


def _msa_col_to_window_idx(col: int, cfg: dict) -> int | None:
    aa_start = cfg["aa_start"]
    aa_end = cfg["aa_end"]
    if col < aa_start or col >= aa_end:
        return None
    return col - aa_start


def _set_band(mask: torch.Tensor, cfg: dict, msa, ref_idx: int, lo: int, hi: int) -> int:
    n = 0
    for aa in range(lo, hi + 1):
        col = _ref_pos_to_msa_col(msa, ref_idx, aa)
        if col is None:
            continue
        wi = _msa_col_to_window_idx(col, cfg)
        if wi is not None and 0 <= wi < len(mask):
            mask[wi] = True
            n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--msa",
        default=str(ROOT / "results" / "bdbv_l_conservation" / "l_msa.fasta"),
    )
    args = ap.parse_args()
    cfg = load_window_config()
    L = int(cfg["max_seq_len"])
    mask = torch.zeros(L, dtype=torch.bool)

    msa_path = Path(args.msa)
    if not msa_path.is_file():
        print(f"Missing MSA {msa_path} — run bdbv_define_l_window.py first")
        sys.exit(1)
    msa = AlignIO.read(msa_path, "fasta")

    ref_idx = 0
    for i, rec in enumerate(msa):
        if REF_EBOV_ACCESSION.split(".")[0] in rec.id or rec.id.startswith("KM034562"):
            ref_idx = i
            break

    meta_bands: list[dict] = []

    def band(name: str, lo: int, hi: int, hot: bool = True):
        if hot:
            n = _set_band(mask, cfg, msa, ref_idx, lo, hi)
        else:
            n = 0
        meta_bands.append({"name": name, "aa_lo": lo, "aa_hi": hi, "hot": hot, "sites_set": n})

    band("D759_hotspot", 754, 764)
    band("RdRp_domain", 540, 850)
    band("N_term_VP35_hub", 1, 380)
    # CTD: last 250 aa of ~2211
    band("CTD_methyltransferase", 2211 - 250, 2211)
    band("GDNQ_core", GDNQ_EXCLUDE[0], GDNQ_EXCLUDE[1], hot=False)
    band("BDBV_2026_K1738", 1738, 1738)
    band("BDBV_2026_Q1770", 1770, 1770)

    # MSA entropy peaks in RdRp band refine
    ent_path = ROOT / "results" / "bdbv_l_conservation" / "conservation.tsv"
    if ent_path.is_file():
        ent = []
        with open(ent_path) as f:
            next(f)
            for line in f:
                _, _, e = line.strip().split("\t")
                ent.append(float(e))
        rd_lo, rd_hi = 540, 850
        seg = [(i, ent[i]) for i in range(len(ent)) if rd_lo <= i + 1 <= rd_hi]
        seg.sort(key=lambda x: -x[1])
        for col, _ in seg[:20]:
            wi = _msa_col_to_window_idx(col, cfg)
            if wi is not None:
                mask[wi] = True

    out_dir = ROOT / "results" / "bdbv_l_mask"
    out_dir.mkdir(parents=True, exist_ok=True)
    pt_path = out_dir / "mut_hotspot_mask_lit.pt"
    torch.save(
        {
            "mut_hotspot_mask": mask.cpu(),
            "meta": {
                "max_seq_len": L,
                "window_config": str(window_config_path()),
                "bands": meta_bands,
                "n_hot": int(mask.sum().item()),
                "gdnq_excluded": list(GDNQ_EXCLUDE),
            },
        },
        pt_path,
    )

    eval_variants = {
        "variants": [
            {"name": "L:K1738N", "species": "bdbv", "aa_ref": 1738, "ref": REF_FJ217161},
            {"name": "L:Q1770R", "species": "bdbv", "aa_ref": 1770, "ref": REF_FJ217161},
            {"name": "L:D759G", "species": "ebov", "aa_ref": 759, "ref": REF_EBOV_ACCESSION},
        ],
        "window_config": cfg,
    }
    (out_dir / "mut_hotspot_mask_eval_variants.json").write_text(
        json.dumps(eval_variants, indent=2) + "\n"
    )
    (out_dir / "mut_hotspot_mask_lit.meta.json").write_text(
        json.dumps({"bands": meta_bands, "n_hot": int(mask.sum().item()), "L": L}, indent=2)
        + "\n"
    )
    print(f"Hotspot mask: {int(mask.sum())}/{L} sites -> {pt_path}")


if __name__ == "__main__":
    main()
