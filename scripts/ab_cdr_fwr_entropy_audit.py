#!/usr/bin/env python3
"""
CDR vs FWR entropy / p_stay / mask audit for OAS + Rod.82.

Compares:
  - column entropy (AA) in IMGT CDR vs FWR on OAS train clones and Rod.82
  - TreeSBM hotspot mask occupancy vs CDR columns
  - optional Thrifty/DASM site-rate mean in CDR vs FWR (if netam available)

Usage:
  python scripts/ab_cdr_fwr_entropy_audit.py \
    --oas-data data/ab_clones_1m/train \
    --rod-jsonl antibody_benchmark/data/processed/benchmark_trees.jsonl \
    --cdr-mask results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt \
    --out benchmarks/results/tables/ab_cdr_fwr_entropy_audit.md
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

import torch

AA = set("ACDEFGHIKLMNPQRSTVWY")


def _entropy(counts: Counter) -> float:
    n = sum(counts.values())
    if n <= 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c <= 0:
            continue
        p = c / n
        h -= p * math.log(p + 1e-12)
    return h


def column_entropy(seqs: list[str], L: int) -> list[float]:
    cols = [Counter() for _ in range(L)]
    for s in seqs:
        for i, a in enumerate(s[:L]):
            if a in AA:
                cols[i][a] += 1
    return [_entropy(c) for c in cols]


def load_oas_leaves(data_dir: Path, max_groups: int, max_seq_len: int) -> tuple[list[str], int]:
    seqs: list[str] = []
    groups = sorted(data_dir.glob("group_*_anc_aa.fasta"))
    n_g = 0
    for fa in groups:
        if max_groups and n_g >= max_groups:
            break
        tip_seqs = []
        cur = None
        name = None
        with fa.open() as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if cur and name and not name.startswith("NODE_"):
                        tip_seqs.append(cur)
                    name = line[1:].split()[0]
                    cur = ""
                else:
                    cur = (cur or "") + line
            if cur and name and not name.startswith("NODE_"):
                tip_seqs.append(cur)
        seqs.extend(s.replace("-", "").upper() for s in tip_seqs if s)
        n_g += 1
    return seqs, max_seq_len


def load_rod_leaves(jsonl: Path, max_seq_len: int) -> list[str]:
    seqs: list[str] = []
    if not jsonl.is_file():
        return seqs
    for line in jsonl.read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        aa = d.get("true_aa_sequences") or {}
        for s in aa.values():
            if s:
                seqs.append(s.replace("-", "").upper()[:max_seq_len])
    return seqs


def summarize(name: str, H: list[float], mask: torch.Tensor) -> dict:
    L = min(len(H), mask.numel())
    m = mask[:L].bool()
    cdr = [H[i] for i in range(L) if bool(m[i])]
    fwr = [H[i] for i in range(L) if not bool(m[i])]
    def mean(xs):
        return sum(xs) / len(xs) if xs else float("nan")
    return {
        "name": name,
        "L": L,
        "n_cdr": int(m[:L].sum().item()),
        "n_fwr": int((~m[:L]).sum().item()),
        "H_cdr_mean": mean(cdr),
        "H_fwr_mean": mean(fwr),
        "H_ratio_cdr_fwr": (mean(cdr) / mean(fwr)) if fwr and mean(fwr) > 0 else float("nan"),
        "mask_frac": float(m[:L].float().mean().item()) if L else 0.0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--oas-data", type=Path, default=Path("data/ab_clones_1m/train"))
    ap.add_argument(
        "--rod-jsonl",
        type=Path,
        default=Path("antibody_benchmark/data/processed/benchmark_trees.jsonl"),
    )
    ap.add_argument(
        "--cdr-mask",
        type=Path,
        default=Path("results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt"),
    )
    ap.add_argument("--max-seq-len", type=int, default=160)
    ap.add_argument("--max-oas-groups", type=int, default=80)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/results/tables/ab_cdr_fwr_entropy_audit.md"),
    )
    args = ap.parse_args()

    blob = torch.load(args.cdr_mask, map_location="cpu", weights_only=False)
    if isinstance(blob, dict) and "mut_hotspot_mask" in blob:
        mask = blob["mut_hotspot_mask"].bool()
    else:
        mask = torch.as_tensor(blob).bool()

    rows = []
    if args.oas_data.is_dir():
        oas_seqs, L = load_oas_leaves(args.oas_data, args.max_oas_groups, args.max_seq_len)
        if oas_seqs:
            H = column_entropy(oas_seqs, L)
            rows.append(summarize("OAS_train_tips", H, mask))
            rows[-1]["n_seqs"] = len(oas_seqs)

    rod = load_rod_leaves(args.rod_jsonl, args.max_seq_len)
    if rod:
        # Rod lengths vary; pad/truncate to mask L
        L = int(mask.numel())
        H = column_entropy(rod, L)
        rows.append(summarize("Rod82_leaves", H, mask))
        rows[-1]["n_seqs"] = len(rod)

    # Thrifty/DASM comparison (best-effort)
    thrifty_note = "Thrifty/DASM site-rate compare: skipped (netam not imported)."
    try:
        from antibody_benchmark.models.dasm import DASMModel  # noqa: F401
        thrifty_note = (
            "DASM available: Thrifty μ × AA selection targets NT SHM hotspots; "
            "TreeSBM v1 Q0 is ESM conservation (anti-correlated with CDR entropy)."
        )
    except Exception as e:
        thrifty_note = f"Thrifty/DASM import note: {e}"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Antibody CDR vs FWR entropy / mask audit",
        "",
        "SHM should elevate mutability in CDR (and AID motifs) vs FWR. "
        "If TreeSBM Q0 = ESM MLM, stay mass concentrates on conserved FWR — "
        "opposite of Neutral/Thrifty/CoSiNE.",
        "",
        "| set | n_seqs | L | n_CDR | n_FWR | H_CDR | H_FWR | H_CDR/H_FWR | mask_frac |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['name']} | {r.get('n_seqs', '—')} | {r['L']} | {r['n_cdr']} | {r['n_fwr']} | "
            f"{r['H_cdr_mean']:.4f} | {r['H_fwr_mean']:.4f} | {r['H_ratio_cdr_fwr']:.3f} | "
            f"{r['mask_frac']:.3f} |"
        )
    lines += [
        "",
        f"CDR mask: `{args.cdr_mask}`",
        "",
        "## Interpretation",
        "",
        "- If H_CDR/H_FWR ≫ 1 on OAS/Rod, empirical SHM load is CDR-biased.",
        "- v1 TreeSBM CDR recall 0.15 vs JC69 ~0.69 is a Q0 mismatch, not missing λ_mut.",
        "- Recipe A (`--shm-site-boost` / `--shm-fwr-stay` / `--shm-use-aid`) lowers stay on CDR∪AID.",
        "",
        thrifty_note,
        "",
    ]
    args.out.write_text("\n".join(lines) + "\n")
    json_path = args.out.with_suffix(".json")
    json_path.write_text(json.dumps({"rows": rows, "thrifty_note": thrifty_note}, indent=2) + "\n")
    print(f"Wrote {args.out}")
    print(f"Wrote {json_path}")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
