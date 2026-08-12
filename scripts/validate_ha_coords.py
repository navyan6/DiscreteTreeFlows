#!/usr/bin/env python3
"""
Validate H3 (or H1) HA columns ↔ paper numbering + write lit hotspot masks.

H3 policy (mirrors COVID spike validation):
  - Paper / WHO / EVEscape use **H3 numbering** (mature HA1; residue 1 after
    the 16-aa signal peptide).
  - Our TreeSBM sequences are **full-length ungapped HA ORF** (signal included,
    L≈566). Same length / pad-to-566 alone does NOT create homology — need
    reference-anchored columns. Empirically group_*_anc_aa are ungapped.
  - Mapping for H3: col_0based = h3_pos + 16 - 1  (orf_1based = h3_pos + 16).
  - Proof: signal motif MKTIIALSYILCLVFA, mature start QKLPGNDNSTATLCL,
    RBS landmarks Y98 / W153 / H183 / Y195, TX/50 F159 / K189 / N225.

Usage:
  python scripts/validate_ha_coords.py \\
      --data data/train --max-seq-len 566 \\
      --out-dir results/flu_mutfreq_vs_lit
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.flu_lit_sites import (  # noqa: E402
    CANONICAL_H1_ORF_LEN,
    H1_ANTIGENIC_SITES_ORF_MET,
    H1_FIELD_GUIDANCE,
    H1_GLOBULAR_HEAD_RANGE,
    H1_MATURE_MOTIF,
    H1_RBS_ORF_MET,
    H1_SIGNAL_LEN_DEFAULT,
    H1_SIGNAL_MOTIF_PREFIX,
    H3_GLOBULAR_HEAD_RANGE,
    H3_SIGNAL_LEN,
    NMICROBIOL201658_H1_PRIMARY,
    NMICROBIOL201658_H3_PRIMARY,
    NMICROBIOL201658_H3_SECONDARY,
    h1_guidance_cols,
    h1_mature_to_col,
    h1_orf_to_col,
    h3_to_col,
    unique_h1_lit_cols,
    unique_h1_primary_positions,
    unique_h3_primary_positions,
)

# TX/50-era landmarks: (label, H3_pos, expected_AA_on_TX50_or_canonical)
# Consensus on modern train MSAs may differ (e.g. 159Y after 2014–15).
H3_LANDMARKS_TX50 = [
    ("signal_end_A16", None, None),  # special-cased
    ("mature_Q1", 1, "Q"),
    ("Y98_RBS", 98, "Y"),
    ("W153_RBS", 153, "W"),
    ("H183_RBS", 183, "H"),
    ("Y195_RBS", 195, "Y"),
    ("F159_TX50", 159, "F"),
    ("K189_TX50", 189, "K"),
    ("N225_TX50", 225, "N"),
]

SIGNAL_MOTIF = "MKTIIALSYILCLVFA"
MATURE_MOTIF = "QKLPGNDNSTATLCL"


def _read_wt(path: Path) -> str:
    if not path.exists():
        return ""
    return "".join(
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith(">")
    )


def _read_fastas(paths: list[Path]) -> list[str]:
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


def _consensus(seqs: list[str], col: int) -> tuple[str, Counter]:
    c: Counter = Counter(s[col] for s in seqs if len(s) > col and s[col] not in "-X*")
    if not c:
        return "-", c
    return c.most_common(1)[0][0], c


def _consensus_seq(seqs: list[str], L: int) -> str:
    cons: list[str] = []
    for i in range(L):
        c = Counter(s[i] for s in seqs if len(s) > i and s[i] not in "-X*")
        cons.append(c.most_common(1)[0][0] if c else "-")
    return "".join(cons)


def validate_h1(args: argparse.Namespace) -> int:
    """Validate H1 pdm09 numbering and write H1-specific lit masks (never overwrite H3)."""
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    signal = args.signal_len if args.signal_len != H3_SIGNAL_LEN else H1_SIGNAL_LEN_DEFAULT
    L = args.max_seq_len

    data = Path(args.data)
    fas = sorted(data.glob("*_anc_aa.fasta"))
    if not fas:
        raise SystemExit(f"No *_anc_aa.fasta under {data}")
    seqs = [s for s in _read_fastas(fas) if len(s) in (L, L - 1, CANONICAL_H1_ORF_LEN)]
    full = [s for s in seqs if len(s) == L]
    if len(full) < 100:
        raise SystemExit(f"Too few L={L} anc seqs under {data}: {len(full)}")

    cons = _consensus_seq(full, L)
    print(f"H1 subtype  data={data}  n_full={len(full)}  signal_len={signal}  L={L}")
    print(f"  cons[:40]={cons[:40]}")

    sig = cons[:signal]
    sig_ok = sig.startswith(H1_SIGNAL_MOTIF_PREFIX[:7]) and len(sig) == signal
    mat = cons[signal : signal + len(H1_MATURE_MOTIF)]
    mat_ok = mat == H1_MATURE_MOTIF
    print(f"  signal[:{signal}]={sig} prefix={H1_SIGNAL_MOTIF_PREFIX} [{'OK' if sig_ok else 'BAD'}]")
    print(f"  mature@{signal}={mat} expect={H1_MATURE_MOTIF} [{'OK' if mat_ok else 'BAD'}]")

    landmarks = []
    # nmicrobiol Sa 153–156 (mature_h1)
    for mut, pos, site, note in NMICROBIOL201658_H1_PRIMARY:
        if pos < 150 or pos > 156:
            continue
        col = h1_mature_to_col(pos, signal)
        aa, counts = _consensus(full, col)
        expect = mut[0]
        ok = aa == expect or (pos == 156 and aa in "NK")  # N156K common
        landmarks.append({
            "mutation": mut,
            "scheme": "mature_h1",
            "pos_1based": pos,
            "orf_1based": pos + signal,
            "col_0based": col,
            "expect_aa": expect,
            "consensus_aa": aa,
            "top_counts": counts.most_common(4),
            "ok": ok,
            "site": site,
            "note": note,
        })
        print(f"  Sa mature#{pos} col={col} expect={expect} got={aa} [{'OK' if ok else 'BAD'}]")

    # Field guidance
    for mut, pos, scheme, note in H1_FIELD_GUIDANCE:
        col = (
            h1_mature_to_col(pos, signal) if scheme == "mature_h1" else h1_orf_to_col(pos)
        )
        aa, counts = _consensus(full, col)
        expect = mut[0]
        alt = mut[-1]
        ok = aa in (expect, alt)  # allow already-drifted consensus
        landmarks.append({
            "mutation": mut,
            "scheme": scheme,
            "pos_1based": pos,
            "col_0based": col,
            "expect_aa": expect,
            "alt_aa": alt,
            "consensus_aa": aa,
            "top_counts": counts.most_common(4),
            "ok": ok,
            "note": note,
        })
        print(
            f"  field {mut} ({scheme}#{pos}) col={col} "
            f"expect∈{{{expect},{alt}}} got={aa} [{'OK' if ok else 'BAD'}]"
        )

    # RBS
    for orf_pos, expect, note in H1_RBS_ORF_MET:
        col = h1_orf_to_col(orf_pos)
        aa, counts = _consensus(full, col)
        ok = aa == expect
        landmarks.append({
            "mutation": f"RBS_{orf_pos}{expect}",
            "scheme": "orf_met",
            "pos_1based": orf_pos,
            "col_0based": col,
            "expect_aa": expect,
            "consensus_aa": aa,
            "top_counts": counts.most_common(3),
            "ok": ok,
            "note": note,
        })
        print(f"  RBS orf#{orf_pos} col={col} expect={expect} got={aa} [{'OK' if ok else 'BAD'}]")

    # Write H1 lit mask (antigenic guidance; NOT full head continuum)
    guidance = h1_guidance_cols(signal)
    lit_mask = torch.zeros(L, dtype=torch.bool)
    lit_cols = []
    lit_labels = []
    for label, col, tag in guidance:
        if 0 <= col < L:
            lit_mask[col] = True
            lit_cols.append(col)
            lit_labels.append({"label": label, "col_0based": col, "tag": tag})
    lit_path = out / "mut_hotspot_mask_h1_nmicrobiol_lit.pt"
    torch.save(
        {
            "mut_hotspot_mask": lit_mask,
            "score": "h1_nmicrobiol_antigenic_guidance",
            "pathogen": "flu_h1n1_ha",
            "signal_len": signal,
            "max_seq_len": L,
            "col_0based": lit_cols,
            "n_sites": int(lit_mask.sum()),
            "labels": lit_labels,
            "nmicrobiol_mature_h1": unique_h1_primary_positions(),
            "antigenic_sites_orf_met": H1_ANTIGENIC_SITES_ORF_MET,
            "field_guidance": [
                {"mutation": m, "pos": p, "scheme": s, "note": n}
                for m, p, s, n in H1_FIELD_GUIDANCE
            ],
            "rbs_orf_met": [
                {"orf_pos": p, "aa": a, "note": n} for p, a, n in H1_RBS_ORF_MET
            ],
            "indexing": (
                f"Dual numbering: mature_h1 col=pos+{signal}-1; "
                f"orf_met col=pos-1. STRICT: not H3 mask, not COVID PMC."
            ),
            "source": (
                "Li et al. Nat Microbiol 2016 (H1 Sa 153–156); "
                "Frontiers 2017 Sa/Sb/Ca/Cb; Matsuzaki/Mountain West field muts"
            ),
        },
        lit_path,
    )
    print(f"Wrote {lit_path} n={int(lit_mask.sum())}")

    # Globular head band (mature_h1)
    lo, hi = H1_GLOBULAR_HEAD_RANGE
    head_mask = torch.zeros(L, dtype=torch.bool)
    for pos in range(lo, hi + 1):
        col = h1_mature_to_col(pos, signal)
        if 0 <= col < L:
            head_mask[col] = True
    head_path = out / "mut_hotspot_mask_h1_globular_head.pt"
    torch.save(
        {
            "mut_hotspot_mask": head_mask,
            "score": "h1_globular_head_band",
            "pathogen": "flu_h1n1_ha",
            "mature_h1_range_1based": [lo, hi],
            "signal_len": signal,
            "max_seq_len": L,
            "note": "H1 globular-head span; separate from H3 head mask and COVID RBD",
        },
        head_path,
    )
    print(f"Wrote {head_path} n={int(head_mask.sum())}")

    union = lit_mask | head_mask
    union_path = out / "mut_hotspot_mask_h1_lit_or_head.pt"
    torch.save(
        {
            "mut_hotspot_mask": union,
            "score": "h1_lit_OR_globular_head",
            "pathogen": "flu_h1n1_ha",
            "n_lit": int(lit_mask.sum()),
            "n_head": int(head_mask.sum()),
            "n_union": int(union.sum()),
            "max_seq_len": L,
            "signal_len": signal,
        },
        union_path,
    )
    print(f"Wrote {union_path} n={int(union.sum())}")

    n_ok = sum(1 for r in landmarks if r.get("ok"))
    proof = {
        "pathogen": "flu_h1n1_ha",
        "indexing": {
            "scheme": "H1 dual: mature_h1 + orf_met vs full-ORF columns",
            "signal_len": signal,
            "mature_to_col_0based": f"col = mature_h1 + {signal} - 1",
            "orf_met_to_col_0based": "col = orf_met - 1",
            "max_seq_len": L,
            "canonical_ha_len": CANONICAL_H1_ORF_LEN,
            "do_not_reuse_h3_signal_16": True,
            "do_not_use_covid_pmc_mask": True,
        },
        "msa": {
            "data": str(data),
            "n_full_len": len(full),
            "consensus_start": cons[:40],
            "signal_ok": sig_ok,
            "mature_ok": mat_ok,
        },
        "landmarks": landmarks,
        "nmicrobiol_mature_153_156": cons[
            h1_mature_to_col(153, signal) : h1_mature_to_col(156, signal) + 1
        ],
        "antigenic_site_consensus_orf_met": {
            name: "".join(cons[h1_orf_to_col(p)] for p in poss)
            for name, poss in H1_ANTIGENIC_SITES_ORF_MET.items()
        },
        "files": {
            "lit_mask": str(lit_path),
            "globular_head_mask": str(head_path),
            "lit_or_head_mask": str(union_path),
        },
        "n_lit_sites": int(lit_mask.sum()),
        "lit_cols_0based": unique_h1_lit_cols(signal),
        "verdict": (
            "PASS"
            if sig_ok and mat_ok and n_ok >= len(landmarks) - 1
            else "FAIL"
        ),
    }
    # Counter not JSON-serializable in landmarks
    for r in proof["landmarks"]:
        if "top_counts" in r:
            r["top_counts"] = [[a, n] for a, n in r["top_counts"]]
    proof_path = out / "COORD_VALIDATION_H1.json"
    proof_path.write_text(json.dumps(proof, indent=2) + "\n")
    print(f"Wrote {proof_path}")
    print("VERDICT:", proof["verdict"], f"({n_ok}/{len(landmarks)} landmarks OK)")
    return 0 if proof["verdict"] == "PASS" else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="data/train",
                    help="Dir with group_*_anc_aa.fasta. Prefer data/h3n2/train or data/h1n1/train on Betty.")
    ap.add_argument("--max-seq-len", type=int, default=566)
    ap.add_argument("--out-dir", default="results/flu_mutfreq_vs_lit")
    ap.add_argument("--wt", default="data/h3n2/wt_ha.txt",
                    help="Full-length HA AA reference (TX/50 AGL07159.1 by default; unused for --subtype h1)")
    ap.add_argument("--full-len", type=int, default=566)
    ap.add_argument("--signal-len", type=int, default=None,
                    help="Default 16 for h3, 17 for h1")
    ap.add_argument("--subtype", choices=["h3", "h1"], default="h3",
                    help="h3=Koel/nmicrobiol H3 mask; h1=H1 antigenic guidance mask")
    args = ap.parse_args()

    if args.signal_len is None:
        args.signal_len = H1_SIGNAL_LEN_DEFAULT if args.subtype == "h1" else H3_SIGNAL_LEN

    if args.subtype == "h1":
        if args.data in ("data/train",):
            args.data = "data/h1n1/train"
        return validate_h1(args)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    signal = args.signal_len

    wt = _read_wt(Path(args.wt))
    if not wt:
        raise SystemExit(f"No WT HA at {args.wt}")

    print(f"WT len={len(wt)}  signal_len={signal}  (expect full_len≈{args.full_len})")
    wt_ok = []
    # Signal peptide
    sig_ok = wt[:signal] == SIGNAL_MOTIF[:signal] if len(wt) >= signal else False
    wt_ok.append({
        "check": "signal_peptide",
        "expect": SIGNAL_MOTIF[:signal],
        "got": wt[:signal],
        "ok": sig_ok,
    })
    print(f"  signal[:{signal}]={wt[:signal]} expect={SIGNAL_MOTIF[:signal]} "
          f"[{'OK' if sig_ok else 'BAD'}]")
    mat_ok = wt[signal : signal + len(MATURE_MOTIF)] == MATURE_MOTIF if len(wt) > signal else False
    wt_ok.append({
        "check": "mature_start_motif",
        "col_0based_start": signal,
        "expect": MATURE_MOTIF,
        "got": wt[signal : signal + len(MATURE_MOTIF)],
        "ok": mat_ok,
    })
    print(f"  mature@{signal}={wt[signal:signal+15]} expect={MATURE_MOTIF} "
          f"[{'OK' if mat_ok else 'BAD'}]")

    for name, h3pos, expect in H3_LANDMARKS_TX50:
        if h3pos is None:
            continue
        col = h3_to_col(h3pos, signal)
        aa = wt[col] if col < len(wt) else "?"
        ok = aa == expect
        wt_ok.append({
            "mutation": name,
            "h3_pos_1based": h3pos,
            "orf_1based": h3pos + signal,
            "col_0based": col,
            "wt_aa": aa,
            "expect": expect,
            "ok": ok,
        })
        print(
            f"  WT {name}: H3#{h3pos} orf={h3pos+signal} col={col} "
            f"aa={aa} expect={expect} [{'OK' if ok else 'BAD'}]"
        )

    fasta_paths = sorted(Path(args.data).glob("group_*_anc_aa.fasta"))
    print(f"anc_aa fastas under {args.data}: {len(fasta_paths)}")
    full_len_proof = None
    len_hist = None
    alignment_audit = None
    if fasta_paths:
        seqs = _read_fastas(fasta_paths)
        len_hist = dict(Counter(len(s) for s in seqs).most_common())
        n_gap_chars = sum(s.count("-") for s in seqs)
        n_seqs_with_gap = sum(1 for s in seqs if "-" in s)
        print(f"n_seqs={len(seqs)}  len_hist_top={list(len_hist.items())[:8]}")
        print(
            f"gap audit: seqs_with_'-'={n_seqs_with_gap}/{len(seqs)}  "
            f"total_gap_chars={n_gap_chars}"
        )
        alignment_audit = {
            "question": (
                "Is pad/truncate to max_seq_len=566 enough for homologous HA columns?"
            ),
            "answer": (
                "No. Same length alone does not define homology. Need reference-"
                "anchored ungapped full-length HA (signal+HA1+HA2) or a gapped MSA. "
                "Empirically anc_aa are left-aligned ungapped; H3# ↔ col only when "
                f"signal_len={signal} and L={args.full_len}."
            ),
            "group_anc_aa_are_gapped_msa": n_gap_chars > 0,
            "n_seqs": len(seqs),
            "n_seqs_with_gap_char": n_seqs_with_gap,
            "total_gap_chars": int(n_gap_chars),
            "max_seq_len_pad_truncate_creates_homology": False,
            "h3_to_col": f"col = h3_pos + {signal} - 1",
        }
        full = [s for s in seqs if len(s) == args.full_len]
        print(f"full-length L={args.full_len}: {len(full)} / {len(seqs)}")
        landmark_msa = []
        if full:
            # Motif locations
            n_sig = sum(1 for s in full[:2000] if s[:signal] == SIGNAL_MOTIF[:signal])
            n_mat = sum(
                1 for s in full[:2000]
                if s[signal : signal + len(MATURE_MOTIF)] == MATURE_MOTIF
            )
            print(f"  motif: signal@{0}={n_sig}/{min(2000,len(full))}  "
                  f"mature@{signal}={n_mat}/{min(2000,len(full))}")
            ids = [
                sum(s[i] == wt[i] for i in range(min(len(s), len(wt)))) / max(len(wt), 1)
                for s in full[:5000]
            ]
            mean_id = sum(ids) / len(ids) if ids else 0.0
            print(f"mean identity to WT (L={args.full_len} sample): {mean_id:.4f}")
            for name, h3pos, expect in H3_LANDMARKS_TX50:
                if h3pos is None:
                    continue
                col = h3_to_col(h3pos, signal)
                cons, counts = _consensus(full, col)
                # 159 may be F (TX/50) or Y/S (post-2014); RBS landmarks should hold
                flexible = name.startswith("F159") or name.startswith("K189") or name.startswith("N225")
                ok = (cons == expect) or (flexible and cons not in ("-", "?"))
                landmark_msa.append({
                    "mutation": name,
                    "h3_pos_1based": h3pos,
                    "col_0based": col,
                    "wt_aa": expect,
                    "consensus_aa": cons,
                    "top_counts": counts.most_common(4),
                    "ok": ok,
                    "note": "polymorphic OK if flexible landmark" if flexible else "",
                })
                print(
                    f"  MSA L={args.full_len} {name}: cons={cons} wt={expect} "
                    f"top={counts.most_common(3)} [{'OK' if ok else 'BAD'}]"
                )
            rbs_ok = all(
                r["ok"] for r in landmark_msa
                if r["mutation"] in ("Y98_RBS", "W153_RBS", "H183_RBS", "Y195_RBS")
            )
            n_check = min(2000, len(full))
            # Signal motif is the hard homology proof. Mature motif can be <90% on
            # mixed data/train (non-H3 groups); require ≥80% or pure-H3 ≥90%.
            signal_ok = n_sig > 0.9 * n_check
            mature_ok = n_mat > 0.8 * n_check
            motif_ok = signal_ok and mature_ok
            verdict = (
                f"PASS: L={args.full_len} columns are full-ORF H3 frame; "
                f"h3_pos = col - {signal} + 1"
                if rbs_ok and motif_ok and mean_id > 0.9
                else "FAIL"
            )
            full_len_proof = {
                "n_full_len": len(full),
                "mean_identity_to_wt": mean_id,
                "n_signal_motif": n_sig,
                "n_mature_motif": n_mat,
                "n_checked": n_check,
                "landmarks": landmark_msa,
                "verdict": verdict,
                "note": (
                    "data/train may mix H3 with other flu groups; prefer "
                    "data/h3n2/train on Betty for a pure-H3 proof. RBS + signal "
                    "landmarks are the decisive checks."
                ),
            }

    # Curated lit masks (store H3 positions + column mask in full-ORF frame)
    lit_h3 = unique_h3_primary_positions()
    lit_mask = torch.zeros(args.max_seq_len, dtype=torch.bool)
    lit_cols = []
    for pos in lit_h3:
        col = h3_to_col(pos, signal)
        if 0 <= col < args.max_seq_len:
            lit_mask[col] = True
            lit_cols.append(col)
    lit_path = out / "mut_hotspot_mask_nmicrobiol_lit.pt"
    torch.save(
        {
            "mut_hotspot_mask": lit_mask,
            "score": "nmicrobiol201658_h3_primary",
            "h3_positions_1based": lit_h3,
            "orf_positions_1based": [p + signal for p in lit_h3],
            "col_0based": lit_cols,
            "signal_len": signal,
            "max_seq_len": args.max_seq_len,
            "indexing": (
                f"H3 numbering (mature); col_0based = h3_pos + {signal} - 1 "
                f"on full-length ungapped HA ORF"
            ),
            "source": "Li et al. Nat Microbiol 2016; nmicrobiol201658 / PMC5087998",
            "mutations": [m for m, _, _, _ in NMICROBIOL201658_H3_PRIMARY],
        },
        lit_path,
    )
    print(f"Wrote {lit_path} n={int(lit_mask.sum())}")

    # Also alias path name used in COVID-style docs
    alias = out / "mut_hotspot_mask_flu_lit.pt"
    torch.save(torch.load(lit_path, map_location="cpu", weights_only=False), alias)
    print(f"Wrote {alias}")

    lo, hi = H3_GLOBULAR_HEAD_RANGE
    head_mask = torch.zeros(args.max_seq_len, dtype=torch.bool)
    for pos in range(lo, hi + 1):
        col = h3_to_col(pos, signal)
        if 0 <= col < args.max_seq_len:
            head_mask[col] = True
    head_path = out / "mut_hotspot_mask_h3_globular_head.pt"
    torch.save(
        {
            "mut_hotspot_mask": head_mask,
            "score": "h3_globular_head_band",
            "h3_range_1based": [lo, hi],
            "signal_len": signal,
            "max_seq_len": args.max_seq_len,
            "note": "TX/50 / CUHK / Kwangju library mutagenesis span (H3 #63–252)",
        },
        head_path,
    )
    print(f"Wrote {head_path} n={int(head_mask.sum())}")

    # Union lit ∪ globular (optional soft region prior)
    union = lit_mask | head_mask
    union_path = out / "mut_hotspot_mask_flu_lit_or_head.pt"
    torch.save(
        {
            "mut_hotspot_mask": union,
            "score": "nmicrobiol_lit_OR_globular_head",
            "n_lit": int(lit_mask.sum()),
            "n_head": int(head_mask.sum()),
            "n_union": int(union.sum()),
            "max_seq_len": args.max_seq_len,
            "signal_len": signal,
        },
        union_path,
    )

    proof = {
        "indexing": {
            "scheme": "H3 numbering (mature HA1) vs full-ORF columns",
            "signal_len": signal,
            "h3_to_col_0based": f"col = h3_pos + {signal} - 1",
            "max_seq_len": args.max_seq_len,
            "canonical_ha_len": args.full_len,
            "F159_col_0based": h3_to_col(159, signal),
            "globular_head_cols_0based": [
                h3_to_col(lo, signal),
                h3_to_col(hi, signal),
            ],
            "pad_truncate_does_not_create_homology": True,
        },
        "alignment_policy": alignment_audit,
        "wt_landmarks": wt_ok,
        "msa_full_length_proof": full_len_proof,
        "len_hist": len_hist,
        "primary_h3_positions": lit_h3,
        "secondary_h3_positions": sorted(
            {p for _, p, _, _ in NMICROBIOL201658_H3_SECONDARY}
        ),
        "files": {
            "lit_mask": str(lit_path),
            "flu_lit_alias": str(alias),
            "globular_head_mask": str(head_path),
            "lit_or_head_mask": str(union_path),
        },
        "h1_transfer_note": (
            "Paper also covers A(H1N1)pdm09 with H1 numbering (Sa 153–156 critical). "
            "Do not apply H3 signal offset to H1; use H1 signal (~17 aa) after "
            "separate landmark validation. H1 primary sites in src/flu_lit_sites.py."
        ),
    }
    proof_path = out / "COORD_VALIDATION.json"
    proof_path.write_text(json.dumps(proof, indent=2))
    print(f"Wrote {proof_path}")
    if full_len_proof:
        print("VERDICT:", full_len_proof["verdict"])
        return 0 if full_len_proof["verdict"].startswith("PASS") else 1
    print("VERDICT: WT landmarks only (no anc_aa under --data)")
    return 0 if all(r.get("ok", True) for r in wt_ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
