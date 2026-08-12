#!/usr/bin/env python3
"""
Compute TRAIN MSA column mutation frequency and compare to Nat Microbiol 2016
HA antigenic sites (nmicrobiol201658).

MSA-select → tree-apply (same design as COVID PMC pipeline):
  - Primary score: per-column frac ≠ consensus on TRAIN MSA (anc_aa / FASTA)
  - Hotspot mask = top-k / top-frac OR curated lit mask from paper sites
  - Literature coords: **H3 numbering**; our columns are full-ORF with 16-aa
    signal → col = h3_pos + 16 - 1

Usage:
  # Local mixed H3 trees (31 groups under data/train):
  python scripts/compute_flu_mut_freq_vs_lit.py \\
      --data data/train --max-seq-len 566 \\
      --out-dir results/flu_mutfreq_vs_lit --hotspot-frac 0.15

  # Betty H3N2 temporal (preferred when trees exist):
  python scripts/compute_flu_mut_freq_vs_lit.py \\
      --data data/h3n2/train --max-seq-len 566 \\
      --out-dir results/flu_mutfreq_vs_lit --hotspot-frac 0.15
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.bridge.site_stats import (  # noqa: E402
    compute_msa_column_mut_freq_bundle,
    compute_tree_mut_freq_bundle,
    hotspot_mask_from_scores,
)
from src.flu_lit_sites import (  # noqa: E402
    H3_GLOBULAR_HEAD_RANGE,
    H3_SIGNAL_LEN,
    NMICROBIOL201658_H1_PRIMARY,
    NMICROBIOL201658_H3_PRIMARY,
    NMICROBIOL201658_H3_SECONDARY,
    h3_to_col,
    unique_h3_primary_positions,
)


def _read_fasta_seqs(paths: list[Path]) -> list[str]:
    seqs: list[str] = []
    for path in paths:
        name = None
        buf: list[str] = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        seqs.append("".join(buf).upper())
                    name = line[1:]
                    buf = []
                else:
                    buf.append(line)
            if name is not None:
                seqs.append("".join(buf).upper())
    return seqs


def _parse_paper_mutations(md_path: Path | None) -> list[tuple[str, int]]:
    if md_path is None or not md_path.exists():
        return []
    text = md_path.read_text(errors="ignore")
    found = sorted(set(re.findall(r"\b([A-Z]\d{2,3}[A-Z])\b", text)))
    out = []
    for m in found:
        pos = int(re.search(r"\d+", m).group())
        out.append((m, pos))
    return out


def _wt(path: Path) -> str:
    if not path.exists():
        return ""
    return "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--max-seq-len", type=int, default=566)
    ap.add_argument("--out-dir", default="results/flu_mutfreq_vs_lit")
    ap.add_argument("--hotspot-topk", type=int, default=None)
    ap.add_argument("--hotspot-frac", type=float, default=None)
    ap.add_argument("--fasta", nargs="*", default=None)
    ap.add_argument("--wt", default="data/h3n2/wt_ha.txt")
    ap.add_argument("--signal-len", type=int, default=H3_SIGNAL_LEN)
    ap.add_argument(
        "--paper-md",
        default=str(
            Path.home()
            / ".cursor/projects/Users-navyanori-Documents-GitHub-DiscreteTreeFlows/uploads/nmicrobiol201658-0.md"
        ),
    )
    ap.add_argument("--skip-tree-diag", action="store_true", default=True)
    ap.add_argument("--with-tree-diag", action="store_true")
    ap.add_argument("--topk-grid", default="9,16,32,64,85,96")
    args = ap.parse_args()

    if args.with_tree_diag:
        args.skip_tree_diag = False
    if args.hotspot_topk is not None and args.hotspot_frac is not None:
        raise SystemExit("Pass only one of --hotspot-topk / --hotspot-frac")
    if args.hotspot_topk is None and args.hotspot_frac is None:
        args.hotspot_frac = 0.15

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    signal = args.signal_len

    data_dir = Path(args.data)
    tree_fasta = sorted(data_dir.glob("group_*_anc_aa.fasta"))
    n_trees = len(tree_fasta)
    fasta_seqs = _read_fasta_seqs([Path(p) for p in args.fasta]) if args.fasta else None

    if n_trees == 0 and not fasta_seqs:
        print(f"ERROR: no group_*_anc_aa.fasta under {args.data} and no --fasta")
        return 2

    print(f"trees(anc_aa)={n_trees}  max_seq_len={args.max_seq_len}  signal={signal}")

    if fasta_seqs is None and tree_fasta:
        print(f"Loading MSA from {len(tree_fasta)} anc_aa FASTAs...")
        msa_seqs = _read_fasta_seqs(tree_fasta)
        print(f"  n_seqs={len(msa_seqs)}")
        msa = compute_msa_column_mut_freq_bundle(None, args.max_seq_len, sequences=msa_seqs)
    elif fasta_seqs and n_trees == 0:
        msa = compute_msa_column_mut_freq_bundle(None, args.max_seq_len, sequences=fasta_seqs)
    else:
        tree_seqs = _read_fasta_seqs(tree_fasta) if tree_fasta else []
        msa = compute_msa_column_mut_freq_bundle(
            None, args.max_seq_len, sequences=tree_seqs + (fasta_seqs or [])
        )

    tree_bundle = None
    if not args.skip_tree_diag and n_trees > 0:
        from src.dataset import TreeDataset

        print("Computing tree-edge diagnostics...")
        dataset = TreeDataset(args.data, max_seq_len=args.max_seq_len)
        tree_bundle = compute_tree_mut_freq_bundle(dataset, args.max_seq_len)

    msa_freq = torch.tensor(msa["msa_mut_freq"], dtype=torch.float32)
    active = msa["n_valid"] > 0
    L_active = int(active.sum())

    if args.hotspot_topk is not None:
        hot = hotspot_mask_from_scores(msa_freq, topk=args.hotspot_topk)
        hot_desc = f"topk={args.hotspot_topk}"
    else:
        hot = hotspot_mask_from_scores(msa_freq, frac=args.hotspot_frac)
        hot_desc = f"frac={args.hotspot_frac}"

    n_hot = int(hot.sum().item())
    print(f"hotspots ({hot_desc}): {n_hot} / {args.max_seq_len} (active {L_active})")

    wt = _wt(Path(args.wt))
    wt_len = len(wt)

    csv_path = out / "msa_mut_freq_by_site.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        header = [
            "col_0based",
            "orf_pos_1based",
            "h3_pos_1based",
            "consensus_aa",
            "wt_aa",
            "msa_mut_freq",
            "n_valid",
            "n_non_consensus",
            "is_hotspot",
        ]
        if tree_bundle is not None:
            header += ["edge_mut_freq", "root_leaf_freq"]
        w.writerow(header)
        for j in range(args.max_seq_len):
            if not active[j]:
                continue
            orf_pos = j + 1
            h3_pos = orf_pos - signal  # may be ≤0 in signal peptide
            wt_aa = wt[j] if j < wt_len else ""
            row = [
                j,
                orf_pos,
                h3_pos if h3_pos >= 1 else "",
                msa["consensus_aa"][j],
                wt_aa,
                f"{msa['msa_mut_freq'][j]:.8f}",
                int(msa["n_valid"][j]),
                int(msa["n_non_consensus"][j]),
                int(bool(hot[j].item())),
            ]
            if tree_bundle is not None:
                row += [
                    f"{tree_bundle['edge_mut_freq'][j]:.8f}",
                    f"{tree_bundle['root_leaf_freq'][j]:.8f}",
                ]
            w.writerow(row)
    print(f"Wrote {csv_path}")

    order = np.argsort(-msa["msa_mut_freq"])
    top_n = min(50, L_active)
    top_path = out / "top_msa_mut_freq_sites.csv"
    with open(top_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "rank", "col_0based", "orf_pos_1based", "h3_pos_1based",
            "consensus_aa", "msa_mut_freq", "is_hotspot",
        ])
        rank = 0
        for j in order:
            if not active[j]:
                continue
            rank += 1
            h3_pos = (j + 1) - signal
            w.writerow([
                rank, int(j), int(j + 1),
                h3_pos if h3_pos >= 1 else "",
                msa["consensus_aa"][j],
                f"{msa['msa_mut_freq'][j]:.8f}",
                int(bool(hot[j].item())),
            ])
            if rank >= top_n:
                break
    print(f"Wrote {top_path}")

    # Literature overlap (H3)
    lit_rows = []
    overlap_primary = []
    lit_h3 = unique_h3_primary_positions()
    lit_set = set(lit_h3)

    def _row(label, mut, h3_pos, domain, note, primary: bool):
        col = h3_to_col(h3_pos, signal)
        in_range = 0 <= col < args.max_seq_len and bool(active[col])
        freq = float(msa["msa_mut_freq"][col]) if in_range else float("nan")
        is_hot = bool(hot[col].item()) if in_range else False
        cons = msa["consensus_aa"][col] if in_range else ""
        wt_aa = wt[col] if col < wt_len else ""
        rec = {
            "set": label,
            "mutation": mut,
            "h3_pos_1based": h3_pos,
            "orf_pos_1based": h3_pos + signal,
            "col_0based": col,
            "domain": domain,
            "note": note,
            "in_alignment": in_range,
            "consensus_aa": cons,
            "wt_aa": wt_aa,
            "msa_mut_freq": freq,
            "in_hotspot_mask": is_hot,
        }
        lit_rows.append(rec)
        if is_hot and primary:
            overlap_primary.append(rec)

    for mut, pos, domain, note in NMICROBIOL201658_H3_PRIMARY:
        _row("primary", mut, pos, domain, note, True)
    for mut, pos, domain, note in NMICROBIOL201658_H3_SECONDARY:
        _row("secondary", mut, pos, domain, note, False)

    rank_of = {}
    r = 0
    for j in order:
        if not active[j]:
            continue
        r += 1
        h3 = (j + 1) - signal
        if h3 >= 1:
            rank_of[h3] = r
    for rec in lit_rows:
        rec["msa_mut_freq_rank"] = rank_of.get(rec["h3_pos_1based"])

    lit_csv = out / "literature_sites_overlap.csv"
    with open(lit_csv, "w", newline="") as f:
        fields = [
            "set", "mutation", "h3_pos_1based", "orf_pos_1based", "col_0based",
            "domain", "note", "in_alignment", "wt_aa", "consensus_aa",
            "msa_mut_freq", "msa_mut_freq_rank", "in_hotspot_mask",
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rec in lit_rows:
            w.writerow(rec)
    print(f"Wrote {lit_csv}")

    paper_parsed = _parse_paper_mutations(Path(args.paper_md))
    (out / "paper_parsed_mutations.json").write_text(
        json.dumps({
            "n": len(paper_parsed),
            "mutations": paper_parsed[:120],
            "note": "Tokens from paywalled upload may be incomplete; curated lists in src/flu_lit_sites.py",
        }, indent=2)
    )

    n_primary = len(lit_set)
    n_overlap = len({r["h3_pos_1based"] for r in overlap_primary})
    print(
        f"Lit overlap (primary ∩ hotspot {hot_desc}): "
        f"{n_overlap}/{n_primary} = {(n_overlap / n_primary if n_primary else 0):.2%}"
    )

    # Precision / recall at top-k (unique primary H3 positions)
    pr_rows = []
    for tok in args.topk_grid.split(","):
        tok = tok.strip()
        if not tok:
            continue
        k = int(tok)
        pred = []
        for j in order:
            if not active[j]:
                continue
            h3 = (j + 1) - signal
            if h3 >= 1:
                pred.append(h3)
            if len(pred) >= k:
                break
        pred_set = set(pred)
        tp = len(pred_set & lit_set)
        precision = tp / len(pred_set) if pred_set else 0.0
        recall = tp / len(lit_set) if lit_set else 0.0
        pr_rows.append({
            "k": k,
            "tp": tp,
            "n_lit": len(lit_set),
            "precision": precision,
            "recall": recall,
            "f1": (
                2 * precision * recall / (precision + recall)
                if (precision + recall) > 0
                else 0.0
            ),
            "hit_positions": sorted(pred_set & lit_set),
            "missed_positions": sorted(lit_set - pred_set),
        })
        print(
            f"  top-{k:3d}: P={precision:.3f}  R={recall:.3f}  "
            f"F1={pr_rows[-1]['f1']:.3f}  hits={sorted(pred_set & lit_set)}"
        )

    pr_path = out / "precision_recall_vs_lit.csv"
    with open(pr_path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "k", "tp", "n_lit", "precision", "recall", "f1",
                "hit_positions", "missed_positions",
            ],
        )
        w.writeheader()
        for row in pr_rows:
            w.writerow({
                **row,
                "hit_positions": ";".join(map(str, row["hit_positions"])),
                "missed_positions": ";".join(map(str, row["missed_positions"])),
            })
    print(f"Wrote {pr_path}")

    # Curated lit mask
    lit_mask = torch.zeros(args.max_seq_len, dtype=torch.bool)
    for pos in lit_h3:
        col = h3_to_col(pos, signal)
        if 0 <= col < args.max_seq_len:
            lit_mask[col] = True
    lit_mask_path = out / "mut_hotspot_mask_nmicrobiol_lit.pt"
    torch.save(
        {
            "mut_hotspot_mask": lit_mask,
            "score": "nmicrobiol201658_h3_primary",
            "h3_positions_1based": lit_h3,
            "orf_positions_1based": [p + signal for p in lit_h3],
            "max_seq_len": args.max_seq_len,
            "signal_len": signal,
            "source": "Li et al. Nat Microbiol 2016; nmicrobiol201658",
            "note": (
                "MSA-select alternative: curated antigenic H3 sites from paper. "
                "Apply via --mut-hotspot-mask (same as COVID PMC lit)."
            ),
        },
        lit_mask_path,
    )
    torch.save(
        torch.load(lit_mask_path, map_location="cpu", weights_only=False),
        out / "mut_hotspot_mask_flu_lit.pt",
    )
    print(f"Wrote {lit_mask_path} (n={int(lit_mask.sum())})")

    lo, hi = H3_GLOBULAR_HEAD_RANGE
    head_mask = torch.zeros(args.max_seq_len, dtype=torch.bool)
    for pos in range(lo, hi + 1):
        col = h3_to_col(pos, signal)
        if 0 <= col < args.max_seq_len and active[col]:
            head_mask[col] = True
    torch.save(
        {
            "mut_hotspot_mask": head_mask,
            "score": "h3_globular_head_band",
            "h3_range_1based": [lo, hi],
            "max_seq_len": args.max_seq_len,
            "signal_len": signal,
        },
        out / "mut_hotspot_mask_h3_globular_head.pt",
    )

    hot_h3 = {(j + 1) - signal for j in range(args.max_seq_len) if hot[j].item()}
    hot_h3 = {p for p in hot_h3 if p >= 1}
    actual_recall = len(hot_h3 & lit_set) / len(lit_set) if lit_set else 0.0
    actual_precision = len(hot_h3 & lit_set) / len(hot_h3) if hot_h3 else 0.0
    recommend = (
        "nmicrobiol_lit_mask" if actual_recall < 0.75 else "msa_mut_freq_topk_or_frac"
    )
    print(
        f"Hotspot vs lit: P={actual_precision:.3f} R={actual_recall:.3f}  "
        f"recommend={recommend}"
    )

    # Persist mut_freq mask
    mask_path = out / "mut_hotspot_mask_mut_freq.pt"
    torch.save(
        {
            "mut_hotspot_mask": hot.cpu(),
            "msa_mut_freq": msa_freq.cpu(),
            "score": "mut_freq",
            "topk": args.hotspot_topk,
            "frac": args.hotspot_frac,
            "max_seq_len": args.max_seq_len,
            "data": args.data,
            "signal_len": signal,
        },
        mask_path,
    )
    print(f"Wrote {mask_path}")

    summary = {
        "design": "MSA-select (msa_mut_freq) → tree-apply (mut_hotspot_mask)",
        "data": args.data,
        "n_trees": n_trees,
        "max_seq_len": args.max_seq_len,
        "active_cols": L_active,
        "hotspot": hot_desc,
        "n_hotspots": n_hot,
        "msa_mut_freq_mean": float(msa_freq[torch.as_tensor(active)].mean()) if L_active else None,
        "msa_mut_freq_max": float(msa_freq.max()),
        "indexing": {
            "scheme": "H3 numbering (mature) vs full-ORF columns",
            "signal_len": signal,
            "h3_to_col": f"col_0based = h3_pos + {signal} - 1",
            "wt_len": wt_len,
            "note": "H3 HA TreeSBM uses max_seq_len=566 (full ORF incl. signal)",
        },
        "literature": {
            "source": "Li et al. Nat Microbiol 2016; nmicrobiol201658 / PMC5087998",
            "n_primary_mutations": len(NMICROBIOL201658_H3_PRIMARY),
            "n_primary_unique_h3_positions": n_primary,
            "n_primary_in_hotspot": n_overlap,
            "primary_in_hotspot": [r["mutation"] for r in overlap_primary],
            "hotspot_precision_vs_lit": actual_precision,
            "hotspot_recall_vs_lit": actual_recall,
            "recommend_mask": recommend,
            "precision_recall_at_k": pr_rows,
            "h1_primary_positions_note": (
                "H1 Sa sites (H1 numbering): "
                + ", ".join(str(p) for _, p, _, _ in NMICROBIOL201658_H1_PRIMARY)
                + " — do not map with H3 signal offset"
            ),
        },
        "files": {
            "by_site": str(csv_path),
            "top_sites": str(top_path),
            "lit_overlap": str(lit_csv),
            "lit_mask": str(lit_mask_path),
            "mut_freq_mask": str(mask_path),
        },
    }
    summary_path = out / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
