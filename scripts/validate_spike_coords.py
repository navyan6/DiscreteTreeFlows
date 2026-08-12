#!/usr/bin/env python3
"""
Validate COVID sequence columns ↔ UniProt P0DTC2 (Wuhan-Hu-1 spike) 1-based indexing.

Policy (length vs homology):
  - Same length alone (pad/truncate to max_seq_len) does NOT define homologous
    spike positions. Column j is only residue j+1 when sequences are
    *column-aligned* to a reference (gapped MSA) or each seq is mapped to P0DTC2.
  - Once MSA-aligned, shared width L makes a spike window = column range; that
    window is 1:1 with P0DTC2 only if the alignment is reference-anchored.
  - Empirically, group_*_anc_aa.fasta are left-aligned ungapped AA strings
    (no '-' gap characters). Landmark checks on L=1273 prove Wuhan-frame
    mapping for full-length seqs; indel-shifted shorter seqs break absolute cols.

Proof target: for full-length ungapped spike (L=1273), column j is residue j+1,
so D614G lives at col 613, N501Y at col 500, RBD ~cols 318–540, etc.

Also (re)writes curated PMC lit / RBD hotspot masks used by train + enrichment.

Usage (Betty preferred; needs group_*_anc_aa.fasta):
  python scripts/validate_spike_coords.py \\
      --data data/covid/train --max-seq-len 1280 \\
      --out-dir results/covid_mutfreq_vs_lit
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

from scripts.compute_msa_mut_freq_vs_lit import (  # noqa: E402
    PMC10142771_PRIMARY,
    PMC10142771_SECONDARY,
)

LANDMARKS = [
    ("D614G", 614, "D"),
    ("N501Y", 501, "N"),
    ("K417N", 417, "K"),
    ("E484K", 484, "E"),
    ("L452R", 452, "L"),
    ("T478K", 478, "T"),
    ("S477N", 477, "S"),
    ("A222V", 222, "A"),
    ("P681R", 681, "P"),
    ("M1", 1, "M"),
]


def _read_wt(path: Path) -> str:
    if not path.exists():
        return ""
    return "".join(
        line.strip() for line in path.read_text().splitlines() if line.strip() and not line.startswith(">")
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="data/covid/train")
    ap.add_argument("--max-seq-len", type=int, default=1280)
    ap.add_argument("--out-dir", default="results/covid_mutfreq_vs_lit")
    ap.add_argument("--wt", default="data/covid/wt.txt")
    ap.add_argument(
        "--evescape",
        default="data/covid/evescape_spike_rbd.pt",
        help="Optional EVEscape blob; reference_seq used as WT fallback / cross-check",
    )
    ap.add_argument("--full-len", type=int, default=1273, help="Canonical Wuhan spike AA length")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    wt = _read_wt(Path(args.wt))
    evescape_ref = ""
    evescape_match = None
    evescape_positions = None
    ev_path = Path(args.evescape)
    if ev_path.exists():
        blob = torch.load(ev_path, map_location="cpu", weights_only=False)
        if isinstance(blob, dict):
            evescape_ref = blob.get("reference_seq") or ""
            evescape_match = blob.get("match_rate")
            evescape_positions = blob.get("positions")
    if not wt and evescape_ref:
        wt = evescape_ref
        print(f"NOTE: using EVEscape reference_seq as WT ({len(wt)} aa)")
    if not wt:
        raise SystemExit(f"No WT spike AA found at {args.wt} and no EVEscape ref")

    print(f"WT len={len(wt)}  (expect {args.full_len})")
    wt_ok = []
    for name, pos, expect in LANDMARKS:
        aa = wt[pos - 1] if pos <= len(wt) else "?"
        ok = aa == expect
        wt_ok.append({"mutation": name, "spike_pos_1based": pos, "col_0based": pos - 1,
                      "wt_aa": aa, "expect": expect, "ok": ok})
        print(f"  WT {name}: pos={pos} col={pos-1} aa={aa} expect={expect} [{'OK' if ok else 'BAD'}]")

    if evescape_ref:
        id_ref = sum(a == b for a, b in zip(wt, evescape_ref)) / max(len(wt), 1)
        print(f"WT vs EVEscape ref identity={id_ref:.6f}  match_rate={evescape_match}")
        if evescape_positions:
            print(f"EVEscape scored positions (1-based): {evescape_positions[0]}..{evescape_positions[-1]} "
                  f"(n={len(evescape_positions)}; RBD band ~331–531)")

    fasta_paths = sorted(Path(args.data).glob("group_*_anc_aa.fasta"))
    print(f"anc_aa fastas under {args.data}: {len(fasta_paths)}")
    full_len_proof = None
    len_hist = None
    indel_note = None
    alignment_audit = None
    if fasta_paths:
        seqs = _read_fastas(fasta_paths)
        len_hist = dict(Counter(len(s) for s in seqs).most_common())
        n_gap_chars = sum(s.count("-") for s in seqs)
        n_seqs_with_gap = sum(1 for s in seqs if "-" in s)
        print(f"n_seqs={len(seqs)}  len_hist_top={list(len_hist.items())[:8]}")
        print(
            f"gap audit: seqs_with_'-'={n_seqs_with_gap}/{len(seqs)}  "
            f"total_gap_chars={n_gap_chars}  "
            f"(0 ⇒ left-aligned ungapped, NOT a gapped MSA)"
        )
        # Sample within-tree length homogeneity
        within_single = within_multi = 0
        for p in fasta_paths[:: max(1, len(fasta_paths) // 40)][:40]:
            gseqs = _read_fastas([p])
            if len({len(s) for s in gseqs}) == 1:
                within_single += 1
            else:
                within_multi += 1
        alignment_audit = {
            "question": (
                "Should sequences be normalized for length first / all same "
                "length to find the spike window?"
            ),
            "answer": (
                "Same length alone (pad/truncate) does NOT define homologous "
                "spike positions. Need a column-aligned MSA (or map each seq "
                "to P0DTC2) so column j is the same residue across sequences/"
                "trees. Once MSA-aligned, shared width L makes a spike window "
                "= column range; 1:1 with P0DTC2 only if reference-anchored. "
                "Landmark AA checks (D614, N501, …) prove mapping for "
                "Wuhan-frame full-length sequences."
            ),
            "group_anc_aa_are_gapped_msa": n_gap_chars > 0,
            "n_seqs": len(seqs),
            "n_seqs_with_gap_char": n_seqs_with_gap,
            "total_gap_chars": int(n_gap_chars),
            "max_seq_len_pad_truncate_creates_homology": False,
            "within_tree_length_homogeneous_sample": {
                "n_groups_sampled": within_single + within_multi,
                "single_length": within_single,
                "multi_length": within_multi,
            },
            "homology_rule": (
                "Absolute column j ↔ P0DTC2 residue j+1 only for sequences "
                "that are Wuhan-frame ungapped full-length (L=1273) or else "
                "explicitly mapped/aligned to that reference. Indel-bearing "
                "shorter left-aligned strings shift all downstream columns."
            ),
        }
        print("ALIGNMENT POLICY:", alignment_audit["answer"])
        full = [s for s in seqs if len(s) == args.full_len]
        print(f"full-length L={args.full_len}: {len(full)} / {len(seqs)}")
        landmark_msa = []
        if full:
            ids = [
                sum(s[i] == wt[i] for i in range(min(len(s), len(wt)))) / len(wt)
                for s in full[:5000]
            ]
            mean_id = sum(ids) / len(ids)
            print(f"mean identity to WT (L={args.full_len} sample): {mean_id:.4f}")
            for name, pos, expect in LANDMARKS:
                cons, counts = _consensus(full, pos - 1)
                # D614G: consensus may be G (fixed) while WT is D — both prove column
                ok = cons == expect or (name == "D614G" and cons in ("D", "G"))
                landmark_msa.append({
                    "mutation": name,
                    "spike_pos_1based": pos,
                    "col_0based": pos - 1,
                    "wt_aa": expect,
                    "consensus_aa_L1273": cons,
                    "top_counts": counts.most_common(4),
                    "ok": ok,
                })
                print(
                    f"  MSA L={args.full_len} {name}: cons={cons} wt={expect} "
                    f"top={counts.most_common(3)} [{'OK' if ok else 'BAD'}]"
                )
            # Motif anchors
            n_g614 = sum(1 for s in full[:2000] if s[611:618] == "YQGVNCT")
            n_d614 = sum(1 for s in full[:2000] if s[611:618] == "YQDVNCT")
            print(f"  motif @ cols 611–617: YQGVNCT={n_g614} YQDVNCT={n_d614} (of ≤2000)")
            full_len_proof = {
                "n_full_len": len(full),
                "mean_identity_to_wt": mean_id,
                "landmarks": landmark_msa,
                "motif_YQGVNCT_at_611": n_g614,
                "motif_YQDVNCT_at_611": n_d614,
                "verdict": (
                    "PASS: L=1273 columns are Wuhan/P0DTC2 frame; "
                    "spike_pos_1based = col_0based + 1"
                    if all(r["ok"] for r in landmark_msa) and mean_id > 0.95
                    else "FAIL"
                ),
            }
        # Indel caveat for shorter left-aligned ungapped seqs
        short = [s for s in seqs if len(s) == 1271]
        if short:
            starts = Counter(s.find("YQGVNCT") for s in short[:500] if "YQGVNCT" in s)
            indel_note = {
                "n_L1271": len(short),
                "YQGVNCT_start_0based": starts.most_common(5),
                "note": (
                    "Shorter sequences are stored ungapped (no '-' for Δ69/70 etc.), "
                    "so absolute column j is NOT homologous to Wuhan residue j+1 after "
                    "upstream deletions. Pad/truncate to max_seq_len=1280 does not fix "
                    "this. Hotspot masks use Wuhan-frame indices (exact on L=1273; "
                    "approximate for indel-shifted trees). Same convention as EVEscape."
                ),
            }
            print("INDEL NOTE:", indel_note["note"])
            print("  L=1271 YQGVNCT starts:", starts.most_common(5), "(expect 611 if unshifted)")

    # Curated masks (Wuhan-frame columns)
    lit_pos = sorted({pos for _, pos, _, _ in PMC10142771_PRIMARY})
    lit_mask = torch.zeros(args.max_seq_len, dtype=torch.bool)
    for pos in lit_pos:
        col = pos - 1
        if 0 <= col < args.max_seq_len:
            lit_mask[col] = True
    lit_path = out / "mut_hotspot_mask_pmc_lit.pt"
    torch.save(
        {
            "mut_hotspot_mask": lit_mask,
            "score": "pmc10142771_primary",
            "spike_positions_1based": lit_pos,
            "max_seq_len": args.max_seq_len,
            "indexing": "spike_pos_1based = col_0based + 1 (Wuhan / P0DTC2 ungapped)",
            "source": "Kumar et al. Viruses 2023; PMC10142771",
            "mutations": [m for m, _, _, _ in PMC10142771_PRIMARY],
        },
        lit_path,
    )
    print(f"Wrote {lit_path} n={int(lit_mask.sum())}")

    rbd_lo, rbd_hi = 319, 541
    rbd_mask = torch.zeros(args.max_seq_len, dtype=torch.bool)
    for pos in range(rbd_lo, rbd_hi + 1):
        col = pos - 1
        if 0 <= col < args.max_seq_len:
            rbd_mask[col] = True
    rbd_path = out / "mut_hotspot_mask_rbd_band.pt"
    torch.save(
        {
            "mut_hotspot_mask": rbd_mask,
            "score": "rbd_band",
            "spike_range_1based": [rbd_lo, rbd_hi],
            "max_seq_len": args.max_seq_len,
            "indexing": "spike_pos_1based = col_0based + 1",
        },
        rbd_path,
    )
    print(f"Wrote {rbd_path} n={int(rbd_mask.sum())}")

    proof = {
        "indexing": {
            "spike_coords": "1-based UniProt P0DTC2 / Wuhan-Hu-1",
            "msa_columns": "0-based; spike_pos = col + 1 for ungapped full-length spike",
            "max_seq_len": args.max_seq_len,
            "canonical_spike_len": args.full_len,
            "D614G_col_0based": 613,
            "RBD_cols_0based": [318, 540],
            "pad_truncate_does_not_create_homology": True,
        },
        "alignment_policy": alignment_audit
        or {
            "question": (
                "Should sequences be normalized for length first / all same "
                "length to find the spike window?"
            ),
            "answer": (
                "Same length alone (pad/truncate) does NOT define homologous "
                "spike positions. Need column-aligned MSA or per-seq map to "
                "P0DTC2; shared L is a window only after alignment, and 1:1 "
                "with P0DTC2 only if reference-anchored."
            ),
            "max_seq_len_pad_truncate_creates_homology": False,
        },
        "wt_landmarks": wt_ok,
        "evescape": {
            "path": str(ev_path) if ev_path.exists() else None,
            "match_rate": evescape_match,
            "ref_identity_to_wt": (
                sum(a == b for a, b in zip(wt, evescape_ref)) / max(len(wt), 1)
                if evescape_ref
                else None
            ),
            "positions_1based_span": (
                [evescape_positions[0], evescape_positions[-1]]
                if evescape_positions
                else None
            ),
        },
        "msa_full_length_proof": full_len_proof,
        "len_hist": len_hist,
        "indel_caveat": indel_note,
        "secondary_lit_positions": sorted({p for _, p, _, _ in PMC10142771_SECONDARY}),
        "files": {"pmc_lit_mask": str(lit_path), "rbd_mask": str(rbd_path)},
    }
    proof_path = out / "COORD_VALIDATION.json"
    proof_path.write_text(json.dumps(proof, indent=2))
    print(f"Wrote {proof_path}")
    if full_len_proof:
        print("VERDICT:", full_len_proof["verdict"])
    return 0 if (not full_len_proof or full_len_proof["verdict"].startswith("PASS")) else 1


if __name__ == "__main__":
    raise SystemExit(main())
