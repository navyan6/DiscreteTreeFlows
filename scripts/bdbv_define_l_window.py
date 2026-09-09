#!/usr/bin/env python3
"""
Pan-filovirus L conservation analysis + fixed training window.

Input: data/bdbv/l_nt/all.fasta
Output:
  results/bdbv_l_conservation/l_msa.fasta
  results/bdbv_l_conservation/conservation.tsv
  results/bdbv_l_conservation/window_config.json
  data/bdbv/l_window/all.fasta  (codon-aligned NT slice)
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scripts.bdbv_common import DATA_ROOT, REF_FJ217161, format_fasta_header, parse_fasta_header, window_config_path

# 2026 BDBV outbreak sites (FJ217161 numbering, 1-based aa)
OUTBREAK_SITES_AA = {"K1738", "Q1770"}
MIN_WIN = 600
MAX_WIN = 900
TARGET_IDENTITY = 0.95


def _shannon_entropy(column: str) -> float:
    counts = Counter(c for c in column if c not in "-.")
    n = sum(counts.values())
    if n == 0:
        return 0.0
    ent = 0.0
    for c in counts.values():
        p = c / n
        ent -= p * math.log2(p)
    return ent


def _identity(column: str) -> float:
    chars = [c for c in column if c not in "-."]
    if not chars:
        return 0.0
    maj = Counter(chars).most_common(1)[0][1]
    return maj / len(chars)


def _map_site_to_column(msa, ref_id: str, aa_1based: int) -> int | None:
    ref = None
    for rec in msa:
        if rec.id.startswith(ref_id.split(".")[0]) or ref_id in rec.id:
            ref = rec
            break
    if ref is None and len(msa) > 0:
        ref = msa[0]
    if ref is None:
        return None
    seq = str(ref.seq)
    pos = 0
    for i, c in enumerate(seq):
        if c == "-":
            continue
        pos += 1
        if pos == aa_1based:
            return i
    return None


def pick_window(n_cols: int, ident: list[float], must_include: set[int]) -> tuple[int, int]:
    best = (0, min(MAX_WIN, n_cols))
    best_score = -1.0
    for start in range(0, max(1, n_cols - MIN_WIN)):
        for end in range(start + MIN_WIN, min(start + MAX_WIN, n_cols) + 1):
            if must_include and not must_include.issubset(set(range(start, end))):
                continue
            seg = ident[start:end]
            if not seg:
                continue
            mean_id = sum(seg) / len(seg)
            length = end - start
            if mean_id >= TARGET_IDENTITY and length >= MIN_WIN:
                score = mean_id * 1000 + length
                if score > best_score:
                    best_score = score
                    best = (start, end)
    if best_score < 0:
        # fallback: include must_include with padding
        if must_include:
            lo = max(0, min(must_include) - 50)
            hi = min(n_cols, max(must_include) + 50)
            hi = max(hi, lo + MIN_WIN)
        else:
            lo, hi = 0, min(MAX_WIN, n_cols)
        return lo, hi
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DATA_ROOT / "l_nt" / "all.fasta"))
    ap.add_argument("--max-seqs", type=int, default=400, help="Subsample for MSA speed")
    args = ap.parse_args()

    inp = Path(args.input)
    if not inp.is_file():
        print(f"Missing {inp} — run bdbv_extract_l.py first")
        sys.exit(1)

    recs = []
    for i, rec in enumerate(SeqIO.parse(inp, "fasta")):
        if i >= args.max_seqs:
            break
        aa = str(Seq(str(rec.seq)).translate())
        recs.append(SeqRecord(Seq(aa), id=rec.id, description=rec.description))
    if len(recs) < 4:
        print(f"Only {len(recs)} sequences — need more data")
        sys.exit(1)

    out_dir = ROOT / "results" / "bdbv_l_conservation"
    out_dir.mkdir(parents=True, exist_ok=True)
    aa_in = out_dir / "_aa_input.fasta"
    SeqIO.write(recs, aa_in, "fasta")
    msa_path = out_dir / "l_msa.fasta"
    with open(msa_path, "w") as msa_out:
        proc = subprocess.run(
            ["mafft", "--auto", "--quiet", str(aa_in)],
            check=False,
            stdout=msa_out,
            stderr=subprocess.PIPE,
            text=True,
        )
    if proc.returncode != 0 or not msa_path.is_file() or msa_path.stat().st_size == 0:
        print(f"MAFFT failed (rc={proc.returncode}): {proc.stderr[:500]}", file=sys.stderr)
        sys.exit(1)
    msa = AlignIO.read(msa_path, "fasta")
    n_cols = msa.get_alignment_length()
    ident = [_identity(msa[:, i]) for i in range(n_cols)]
    ent = [_shannon_entropy(msa[:, i]) for i in range(n_cols)]

    tsv = out_dir / "conservation.tsv"
    with open(tsv, "w") as f:
        f.write("col\tidentity\tentropy\n")
        for i in range(n_cols):
            f.write(f"{i}\t{ident[i]:.4f}\t{ent[i]:.4f}\n")

    must_cols: set[int] = set()
    for site_name, aa_pos in (("K1738", 1738), ("Q1770", 1770)):
        col = _map_site_to_column(msa, REF_FJ217161, aa_pos)
        if col is not None:
            must_cols.add(col)

    aa_start, aa_end = pick_window(n_cols, ident, must_cols)
    aa_len = aa_end - aa_start
    nt_start = aa_start * 3
    nt_end = aa_end * 3

    cfg = {
        "max_seq_len": aa_len,
        "aa_start": aa_start,
        "aa_end": aa_end,
        "nt_start": nt_start,
        "nt_end": nt_end,
        "target_identity": TARGET_IDENTITY,
        "must_include_cols": sorted(must_cols),
        "mean_identity_window": sum(ident[aa_start:aa_end]) / max(1, aa_len),
    }
    window_config_path().parent.mkdir(parents=True, exist_ok=True)
    window_config_path().write_text(json.dumps(cfg, indent=2) + "\n")
    print(f"Window aa [{aa_start}:{aa_end}] len={aa_len} mean_id={cfg['mean_identity_window']:.3f}")

    win_dir = DATA_ROOT / "l_window"
    win_dir.mkdir(parents=True, exist_ok=True)
    win_path = win_dir / "all.fasta"
    n_out = 0
    with open(win_path, "w") as out_f:
        for rec in SeqIO.parse(inp, "fasta"):
            acc, date, country, species, source = parse_fasta_header(rec)
            nt = str(rec.seq)
            if len(nt) < nt_end:
                continue
            slice_nt = nt[nt_start:nt_end]
            if len(slice_nt) % 3 != 0:
                continue
            out_f.write(format_fasta_header(acc, date, country, species, source) + f"\n{slice_nt}\n")
            n_out += 1
    print(f"Wrote {n_out} windowed L CDS -> {win_path}")


if __name__ == "__main__":
    main()
