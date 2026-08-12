#!/usr/bin/env python3
"""
Compute TRAIN MSA column mutation frequency and compare to PMC10142771 spike sites.

MSA-select → tree-apply design:
  - Primary score: per-column frac of sequences ≠ consensus/modal AA on the TRAIN MSA
    (from TreeDataset node sequences, or optional AA FASTA).
  - Hotspot mask = top-k / top-frac of that score (same API as train --mut-hotspot-*).
  - Secondary diagnostics: parent→child edge mut-freq and root→leaf mut-freq on trees.
  - Literature: UniProt P0DTC2 / Wuhan spike 1-based residue indices from Kumar et al.
    Viruses 2023 (PMC10142771). Our MSA columns are 0-based; for ungapped spike AA
    aligned to Wuhan, spike_pos = column_index + 1 (pad columns beyond 1273 ignored).

Usage:
  # COVID trees on Betty (preferred):
  python scripts/compute_msa_mut_freq_vs_lit.py \\
      --data data/covid/train --max-seq-len 1280 \\
      --out-dir results/covid_mutfreq_vs_lit --hotspot-frac 0.15

  # Flu smoke (local):
  python scripts/compute_msa_mut_freq_vs_lit.py \\
      --data data/train --max-seq-len 566 --out-dir results/flu_mutfreq_smoke \\
      --skip-lit
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
from src.dataset import TreeDataset  # noqa: E402

# Primary mutations highlighted in PMC10142771 abstract + dedicated sections.
# Coordinates: 1-based Wuhan-Hu-1 / UniProt P0DTC2 spike residue index.
PMC10142771_PRIMARY = [
    # Abstract / study-selected set (+ D614G reference)
    ("A222V", 222, "NTD", "study-selected"),
    ("N439K", 439, "RBD", "study-selected"),
    ("N501Y", 501, "RBD/RBM", "study-selected / Alpha"),
    ("L452R", 452, "RBD", "study-selected / Delta"),
    ("Y453F", 453, "RBD", "study-selected / Alpha"),
    ("E484K", 484, "RBD/RBM", "study-selected / Beta"),
    ("K417N", 417, "RBD", "study-selected / Beta"),
    ("T478K", 478, "RBD", "study-selected / Delta"),
    ("L981F", 981, "S2", "study-selected"),
    ("L212I", 212, "NTD", "study-selected"),
    ("N856K", 856, "S2", "study-selected"),
    ("T547K", 547, "S1", "study-selected"),
    ("G496S", 496, "RBD", "study-selected"),
    ("Y369C", 369, "RBD", "study-selected"),
    ("D614G", 614, "S1", "global reference"),
    ("S477N", 477, "RBD/RBM", "section 3.3"),
]

# Additional VoC / RBD sites discussed in the paper body (secondary lit set).
PMC10142771_SECONDARY = [
    ("K417T", 417, "RBD", "Gamma"),
    ("E484Q", 484, "RBD", "Delta Plus context"),
    ("E484A", 484, "RBD", "mentioned"),
    ("P681R", 681, "furin loop", "Delta context"),
    ("G339D", 339, "RBD", "mentioned"),
    ("G446S", 446, "RBD", "mentioned"),
    ("N440K", 440, "RBD", "mentioned"),
    ("Q493R", 493, "RBD", "mentioned"),
    ("Q498R", 498, "RBD", "mentioned"),
    ("Y505H", 505, "RBD", "mentioned"),
    ("S371L", 371, "RBD", "mentioned"),
    ("S373P", 373, "RBD", "mentioned"),
    ("S375F", 375, "RBD", "mentioned"),
]


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
    """Extract A222V-style tokens from the uploaded PMC markdown (sanity cross-check)."""
    if md_path is None or not md_path.exists():
        return []
    text = md_path.read_text(errors="ignore")
    found = sorted(set(re.findall(r"\b([A-Z]\d{2,4}[A-Z])\b", text)))
    out = []
    for m in found:
        pos = int(re.search(r"\d+", m).group())
        out.append((m, pos))
    return out


def _spike_wt(path: Path) -> str:
    if not path.exists():
        return ""
    return "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="TreeDataset dir (group_*_rooted.nwk + _anc_aa.fasta)")
    ap.add_argument("--max-seq-len", type=int, default=1280)
    ap.add_argument("--out-dir", default="results/covid_mutfreq_vs_lit")
    ap.add_argument("--hotspot-topk", type=int, default=None)
    ap.add_argument(
        "--hotspot-frac",
        type=float,
        default=None,
        help="Top fraction of columns by msa_mut_freq (default 0.15 if neither topk nor frac set)",
    )
    ap.add_argument(
        "--fasta",
        nargs="*",
        default=None,
        help="Optional AA FASTA(s) instead of / in addition to tree node seqs for MSA score",
    )
    ap.add_argument("--wt", default="data/covid/wt.txt", help="Wuhan spike AA (P0DTC2) for index check")
    ap.add_argument(
        "--paper-md",
        default=str(
            Path.home()
            / ".cursor/projects/Users-navyanori-Documents-GitHub-DiscreteTreeFlows/uploads/PMC10142771-0.md"
        ),
    )
    ap.add_argument("--skip-lit", action="store_true", help="Skip PMC overlap (e.g. flu smoke)")
    ap.add_argument(
        "--skip-tree-diag",
        action="store_true",
        default=True,
        help="Skip secondary tree-edge diagnostics (default on for speed)",
    )
    ap.add_argument(
        "--with-tree-diag",
        action="store_true",
        help="Compute parent→child / root→leaf mut-freq diagnostics (slow)",
    )
    ap.add_argument(
        "--prefer-lit-mask",
        action="store_true",
        help="Also write curated hotspot mask from unique PMC primary spike positions "
             "(preferred when MSA top-frac diverges from antigenic lit sites).",
    )
    ap.add_argument(
        "--topk-grid",
        default="16,32,64,96,192",
        help="Comma-separated top-k values for precision/recall vs lit primary positions",
    )
    args = ap.parse_args()

    if args.with_tree_diag:
        args.skip_tree_diag = False

    if args.hotspot_topk is not None and args.hotspot_frac is not None:
        raise SystemExit("Pass only one of --hotspot-topk / --hotspot-frac")
    if args.hotspot_topk is None and args.hotspot_frac is None:
        args.hotspot_frac = 0.15

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    data_dir = Path(args.data)
    tree_fasta = sorted(data_dir.glob("group_*_anc_aa.fasta"))
    n_trees = len(tree_fasta)
    fasta_seqs = _read_fasta_seqs([Path(p) for p in args.fasta]) if args.fasta else None

    if n_trees == 0 and not fasta_seqs:
        print(f"ERROR: no group_*_anc_aa.fasta under {args.data} and no --fasta provided")
        return 2

    print(f"trees(anc_aa)={n_trees}  max_seq_len={args.max_seq_len}  fasta_seqs={0 if not fasta_seqs else len(fasta_seqs)}")

    # Fast primary path: read anc_aa FASTAs only (MSA-select; no TreeDataset load).
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
        print("Computing secondary tree-edge / root→leaf diagnostics via TreeDataset...")
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
    print(f"hotspots ({hot_desc}): {n_hot} / {args.max_seq_len} (active cols {L_active})")

    wt = _spike_wt(Path(args.wt))
    wt_len = len(wt)

    # Per-site CSV
    csv_path = out / "msa_mut_freq_by_site.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        header = [
            "col_0based",
            "spike_pos_1based",
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
            spike_pos = j + 1
            wt_aa = wt[j] if j < wt_len else ""
            row = [
                j,
                spike_pos,
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

    # Top sites
    order = np.argsort(-msa["msa_mut_freq"])
    top_n = min(50, L_active)
    top_path = out / "top_msa_mut_freq_sites.csv"
    with open(top_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "col_0based", "spike_pos_1based", "consensus_aa", "msa_mut_freq", "is_hotspot"])
        rank = 0
        for j in order:
            if not active[j]:
                continue
            rank += 1
            w.writerow(
                [
                    rank,
                    int(j),
                    int(j + 1),
                    msa["consensus_aa"][j],
                    f"{msa['msa_mut_freq'][j]:.8f}",
                    int(bool(hot[j].item())),
                ]
            )
            if rank >= top_n:
                break
    print(f"Wrote {top_path}")

    # Literature overlap
    lit_rows = []
    overlap_primary = []
    overlap_secondary = []
    if not args.skip_lit:
        hot_spike_pos = {j + 1 for j in range(args.max_seq_len) if hot[j].item() and active[j]}

        def _row(label, mut, pos, domain, note, primary: bool):
            col = pos - 1
            in_range = 0 <= col < args.max_seq_len and bool(active[col])
            freq = float(msa["msa_mut_freq"][col]) if in_range else float("nan")
            is_hot = bool(hot[col].item()) if in_range else False
            cons = msa["consensus_aa"][col] if in_range else ""
            wt_aa = wt[col] if col < wt_len else ""
            edge = (
                float(tree_bundle["edge_mut_freq"][col])
                if tree_bundle is not None and in_range
                else float("nan")
            )
            rec = {
                "set": label,
                "mutation": mut,
                "spike_pos_1based": pos,
                "col_0based": col,
                "domain": domain,
                "note": note,
                "in_alignment": in_range,
                "consensus_aa": cons,
                "wt_aa": wt_aa,
                "msa_mut_freq": freq,
                "edge_mut_freq": edge,
                "in_hotspot_mask": is_hot,
            }
            lit_rows.append(rec)
            if is_hot:
                (overlap_primary if primary else overlap_secondary).append(rec)

        for mut, pos, domain, note in PMC10142771_PRIMARY:
            _row("primary", mut, pos, domain, note, True)
        for mut, pos, domain, note in PMC10142771_SECONDARY:
            _row("secondary", mut, pos, domain, note, False)

        # ranks of primary sites
        rank_of = {}
        r = 0
        for j in order:
            if not active[j]:
                continue
            r += 1
            rank_of[j + 1] = r
        for rec in lit_rows:
            rec["msa_mut_freq_rank"] = rank_of.get(rec["spike_pos_1based"])

        lit_csv = out / "literature_sites_overlap.csv"
        with open(lit_csv, "w", newline="") as f:
            fields = [
                "set",
                "mutation",
                "spike_pos_1based",
                "col_0based",
                "domain",
                "note",
                "in_alignment",
                "wt_aa",
                "consensus_aa",
                "msa_mut_freq",
                "msa_mut_freq_rank",
                "edge_mut_freq",
                "in_hotspot_mask",
            ]
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for rec in lit_rows:
                w.writerow(rec)
        print(f"Wrote {lit_csv}")

        paper_parsed = _parse_paper_mutations(Path(args.paper_md))
        (out / "paper_parsed_mutations.json").write_text(
            json.dumps({"n": len(paper_parsed), "mutations": paper_parsed[:80]}, indent=2)
        )

        n_primary = len({r["spike_pos_1based"] for r in lit_rows if r["set"] == "primary"})
        n_overlap = len({r["spike_pos_1based"] for r in overlap_primary})
        print(
            f"Lit overlap (primary ∩ hotspot {hot_desc}): "
            f"{n_overlap}/{n_primary} unique positions = "
            f"{(n_overlap / n_primary if n_primary else 0):.2%}"
        )
        for rec in sorted(overlap_primary, key=lambda r: -r["msa_mut_freq"]):
            print(
                f"  {rec['mutation']:6s} spike={rec['spike_pos_1based']:4d}  "
                f"msa_freq={rec['msa_mut_freq']:.4f}  rank={rec['msa_mut_freq_rank']}"
            )

        # Precision / recall at top-k (unique primary spike positions as positives)
        lit_pos = sorted({pos for _, pos, _, _ in PMC10142771_PRIMARY})
        lit_set = set(lit_pos)
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
                pred.append(j + 1)
                if len(pred) >= k:
                    break
            pred_set = set(pred)
            tp = len(pred_set & lit_set)
            precision = tp / len(pred_set) if pred_set else 0.0
            recall = tp / len(lit_set) if lit_set else 0.0
            pr_rows.append(
                {
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
                }
            )
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
                w.writerow(
                    {
                        **row,
                        "hit_positions": ";".join(map(str, row["hit_positions"])),
                        "missed_positions": ";".join(map(str, row["missed_positions"])),
                    }
                )
        print(f"Wrote {pr_path}")

        # Curated literature mask (unique primary positions) — preferred train option
        # when raw MSA freq alone under-covers antigenic RBD sites.
        lit_mask = torch.zeros(args.max_seq_len, dtype=torch.bool)
        for pos in lit_pos:
            col = pos - 1
            if 0 <= col < args.max_seq_len:
                lit_mask[col] = True
        lit_mask_path = out / "mut_hotspot_mask_pmc_lit.pt"
        torch.save(
            {
                "mut_hotspot_mask": lit_mask,
                "score": "pmc10142771_primary",
                "spike_positions_1based": lit_pos,
                "max_seq_len": args.max_seq_len,
                "source": "Kumar et al. Viruses 2023; PMC10142771",
                "note": (
                    "MSA-select alternative: curated antigenic/RBD sites from paper. "
                    "Apply on tree bridge via mut_hotspot_mask (same as mut_freq/entropy)."
                ),
            },
            lit_mask_path,
        )
        print(f"Wrote {lit_mask_path} (n={int(lit_mask.sum())})")

        # RBD column band from paper (residues 333–528 cited near A222V; RBD ~319–541 standard)
        rbd_lo, rbd_hi = 319, 541  # inclusive 1-based; common Wuhan RBD span
        rbd_mask = torch.zeros(args.max_seq_len, dtype=torch.bool)
        for pos in range(rbd_lo, rbd_hi + 1):
            col = pos - 1
            if 0 <= col < args.max_seq_len and active[col]:
                rbd_mask[col] = True
        torch.save(
            {
                "mut_hotspot_mask": rbd_mask,
                "score": "rbd_band",
                "spike_range_1based": [rbd_lo, rbd_hi],
                "max_seq_len": args.max_seq_len,
            },
            out / "mut_hotspot_mask_rbd_band.pt",
        )

        # Recommend: if frac hotspot recall < 0.75 on primary lit, prefer lit mask
        frac_recall = next(
            (r["recall"] for r in pr_rows if r["k"] == n_hot),
            (n_overlap / n_primary if n_primary else 0.0),
        )
        # Also evaluate the actual hotspot size against lit
        hot_pos = {j + 1 for j in range(args.max_seq_len) if hot[j].item()}
        actual_recall = len(hot_pos & lit_set) / len(lit_set) if lit_set else 0.0
        actual_precision = len(hot_pos & lit_set) / len(hot_pos) if hot_pos else 0.0
        recommend = (
            "pmc_lit_mask"
            if actual_recall < 0.75
            else "msa_mut_freq_topk_or_frac"
        )
        print(
            f"Hotspot vs lit: P={actual_precision:.3f} R={actual_recall:.3f}  "
            f"recommend={recommend}"
        )

    summary = {
        "design": "MSA-select (msa_mut_freq = frac ≠ consensus) → tree-apply (mut_hotspot_mask)",
        "data": args.data,
        "n_trees": n_trees,
        "max_seq_len": args.max_seq_len,
        "active_cols": L_active,
        "hotspot": hot_desc,
        "n_hotspots": n_hot,
        "msa_mut_freq_mean": float(msa_freq[torch.as_tensor(active)].mean()) if L_active else None,
        "msa_mut_freq_max": float(msa_freq.max()),
        "indexing": {
            "spike_coords": "1-based UniProt P0DTC2 / Wuhan-Hu-1",
            "msa_columns": "0-based; spike_pos = col + 1 when AA MSA is ungapped Wuhan-frame",
            "wt_len": wt_len,
            "note": "COVID TreeSBM uses max_seq_len=1280 (pads past spike 1273)",
        },
        "literature": {
            "source": "Kumar et al. Viruses 2023; PMC10142771",
            "n_primary_mutations": len(PMC10142771_PRIMARY),
            "n_primary_unique_positions": len({pos for _, pos, _, _ in PMC10142771_PRIMARY}),
            "n_primary_in_hotspot": len({r["spike_pos_1based"] for r in overlap_primary}),
            "primary_in_hotspot": [r["mutation"] for r in overlap_primary],
            "n_secondary_in_hotspot": len(overlap_secondary),
            "secondary_in_hotspot": [r["mutation"] for r in overlap_secondary],
            "hotspot_precision_vs_lit": actual_precision,
            "hotspot_recall_vs_lit": actual_recall,
            "recommend_mask": recommend,
            "precision_recall_at_k": pr_rows,
        }
        if not args.skip_lit
        else None,
        "files": {
            "by_site": str(csv_path),
            "top_sites": str(top_path),
        },
    }
    summary_path = out / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {summary_path}")

    # Persist hotspot mask for train reuse / inspection
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
        },
        mask_path,
    )
    print(f"Wrote {mask_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
