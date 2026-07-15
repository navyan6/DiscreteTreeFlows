#!/usr/bin/env python3
"""
Convert an EVEscape per-substitution CSV into a static [max_seq_len, 20] score
tensor aligned to the TreeSBM H3N2-HA column frame.

EVEscape is used as an EXTERNAL EVALUATOR of generated lineages (score escape of
generated mutations), NOT baked into the reference process R0 — the Schrodinger
bridge reference stays pure ESM. This script produces the aligned per-(pos, aa)
score table the evaluator looks up; for scoring, prefer raw scores
(--no-standardize).

Why alignment (not `position - 1`):
    EVEscape scores use "H3 numbering" (mature HA1, residue 1 = first residue
    after the 16-aa signal peptide), while our sequences are full-length HA
    (signal peptide included, columns 0..565). H3 numbering also carries
    structural insertions (e.g. 133a) that break linear indexing. So we
    reconstruct EVEscape's wild-type sequence from the CSV and align it to our
    reference sequence to derive `evescape_position -> alignment_column`
    empirically, which is convention-agnostic and self-checking (we report the
    WT match rate at mapped columns).

Output: {"scores": Tensor[L,20], "reference_seq": str, "positions": list[int],
         "match_rate": float, "standardized": bool, "mean": float, "std": float}

The column frame is shared across all master_h3n2 groups (HA is indel-free;
verified group roots are >98% identical column-for-column), so a single tensor
applies to every H3N2 HA group.
"""

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import torch
from Bio import SeqIO

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_VOCAB)}

# Flexible column-name detection for varied EVEscape CSV exports.
POS_ALIASES   = ["position", "pos", "i", "site", "wt_pos"]
WT_ALIASES    = ["wildtype_aa", "wildtype", "wt_aa", "wt", "wildtype_res", "aa_wt"]
MUT_ALIASES   = ["mutant_aa", "mutant", "mut_aa", "mut", "mutation_aa", "aa_mut"]
SCORE_ALIASES = ["evescape_score", "evescape", "score", "escape_score", "escape"]


def _pick_col(header: list[str], aliases: list[str], override: str | None) -> str:
    if override:
        if override not in header:
            raise SystemExit(f"Column '{override}' not in CSV header {header}")
        return override
    lower = {h.lower(): h for h in header}
    for a in aliases:
        if a in lower:
            return lower[a]
    raise SystemExit(
        f"Could not auto-detect a column among {aliases} in header {header}. "
        f"Pass an explicit --*-col."
    )


def load_reference_seq(args) -> str:
    """Full-length HA reference defining the 566-column frame."""
    if args.ref_seq:
        return args.ref_seq.strip().upper()

    from src.dataset import parse_newick

    g = args.ref_from_group
    anc = ROOT / args.data / f"group_{g:03d}_anc_aa.fasta"
    seqs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(anc, "fasta")}
    # Prefer the tree root; fall back to any full-length record.
    try:
        root_id, _, _, _ = parse_newick(str(ROOT / args.data / f"group_{g:03d}_rooted.nwk"))
        ref = seqs.get(root_id, "")
    except Exception:
        ref = ""
    if len(ref) != args.max_seq_len:
        full = [s for s in seqs.values() if len(s) == args.max_seq_len]
        if not full:
            raise SystemExit(
                f"No length-{args.max_seq_len} sequence found in {anc}; "
                f"pass --ref-seq explicitly."
            )
        ref = full[0]
    return ref


def read_csv_rows(args):
    with open(args.csv_path, newline="") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames or []
        pos_c   = _pick_col(header, POS_ALIASES,   args.pos_col)
        wt_c    = _pick_col(header, WT_ALIASES,    args.wt_col)
        mut_c   = _pick_col(header, MUT_ALIASES,   args.mut_col)
        score_c = _pick_col(header, SCORE_ALIASES, args.score_col)
        print(f"Columns: pos={pos_c} wt={wt_c} mut={mut_c} score={score_c}")
        rows = []
        for r in reader:
            try:
                pos = int(float(r[pos_c]))
                wt  = r[wt_c].strip().upper()
                mut = r[mut_c].strip().upper()
                sc  = float(r[score_c])
            except (ValueError, KeyError):
                continue
            if len(wt) == 1 and len(mut) == 1:
                rows.append((pos, wt, mut, sc))
    return rows


def build_evescape_wt(rows) -> tuple[str, dict[int, int]]:
    """Reconstruct the EVEscape WT sequence and a position->wt_string_index map."""
    wt_by_pos: dict[int, str] = {}
    for pos, wt, _mut, _sc in rows:
        if wt in AA_TO_IDX:
            wt_by_pos.setdefault(pos, wt)
    positions = sorted(wt_by_pos)
    wt_seq = "".join(wt_by_pos[p] for p in positions)
    pos_to_wtidx = {p: i for i, p in enumerate(positions)}
    return wt_seq, pos_to_wtidx


def align_positions(evescape_wt: str, reference: str) -> dict[int, int]:
    """
    Map EVEscape-WT string index -> reference column via global protein alignment.
    Returns {wt_string_index: reference_column}.
    """
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    aln = aligner.align(evescape_wt, reference)[0]

    mapping: dict[int, int] = {}
    # aln.aligned: pairs of (start,end) blocks in seqA (evescape_wt) and seqB (reference)
    blocks_a, blocks_b = aln.aligned
    for (a0, a1), (b0, b1) in zip(blocks_a, blocks_b):
        for k in range(a1 - a0):
            mapping[a0 + k] = b0 + k
    return mapping


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv-path", required=True, help="EVEscape per-substitution CSV")
    p.add_argument("--output", default="data/evescape_h3n2_ha.pt")
    p.add_argument("--max-seq-len", type=int, default=566)
    # reference frame
    p.add_argument("--data", default="data/train")
    p.add_argument("--ref-from-group", type=int, default=1,
                   help="master_h3n2 group whose root defines the 566-col frame")
    p.add_argument("--ref-seq", default=None,
                   help="explicit full-length reference sequence (overrides --ref-from-group)")
    # column overrides (else auto-detected)
    p.add_argument("--pos-col", default=None)
    p.add_argument("--wt-col", default=None)
    p.add_argument("--mut-col", default=None)
    p.add_argument("--score-col", default=None)
    # normalization
    p.add_argument("--standardize", action="store_true", default=True,
                   help="z-score non-zero entries so magnitudes match logit scale")
    p.add_argument("--no-standardize", dest="standardize", action="store_false")
    p.add_argument("--min-match-rate", type=float, default=0.90,
                   help="abort if WT match rate at mapped columns is below this")
    args = p.parse_args()

    L = args.max_seq_len
    reference = load_reference_seq(args)
    print(f"Reference length: {len(reference)}  head: {reference[:24]}")

    rows = read_csv_rows(args)
    if not rows:
        raise SystemExit("No usable rows parsed from CSV.")
    print(f"Parsed {len(rows)} (pos,wt,mut,score) rows")

    evescape_wt, pos_to_wtidx = build_evescape_wt(rows)
    print(f"EVEscape WT reconstructed: {len(evescape_wt)} residues")

    wtidx_to_col = align_positions(evescape_wt, reference)

    # Diagnostics: does the mapped column actually carry the EVEscape WT residue?
    matched = total = 0
    for pos, widx in pos_to_wtidx.items():
        col = wtidx_to_col.get(widx)
        if col is None:
            continue
        total += 1
        if col < len(reference) and reference[col] == evescape_wt[widx]:
            matched += 1
    match_rate = matched / total if total else 0.0
    print(f"WT match rate at mapped columns: {match_rate:.3f} ({matched}/{total})")
    if match_rate < args.min_match_rate:
        raise SystemExit(
            f"Match rate {match_rate:.3f} < {args.min_match_rate}. The EVEscape "
            f"reference likely does not correspond to the master_h3n2 frame — "
            f"check --ref-seq / EVEscape reference before proceeding."
        )

    scores = torch.zeros(L, 20, dtype=torch.float32)
    filled = skipped = 0
    for pos, wt, mut, sc in rows:
        if mut not in AA_TO_IDX or wt == mut:
            skipped += 1
            continue
        widx = pos_to_wtidx.get(pos)
        col = wtidx_to_col.get(widx) if widx is not None else None
        if col is None or col >= L:
            skipped += 1
            continue
        scores[col, AA_TO_IDX[mut]] = sc
        filled += 1
    print(f"Filled {filled} (col,aa) entries; skipped {skipped}")

    mean = std = 0.0
    if args.standardize:
        nz = scores != 0
        if nz.any():
            vals = scores[nz]
            mean = float(vals.mean())
            std = float(vals.std()) or 1.0
            scores[nz] = (vals - mean) / std
            print(f"Standardized non-zero entries: mean={mean:.4f} std={std:.4f}")

    out = ROOT / args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "scores": scores,
            "reference_seq": reference,
            "positions": sorted(pos_to_wtidx),
            "match_rate": match_rate,
            "standardized": bool(args.standardize),
            "mean": mean,
            "std": std,
        },
        out,
    )
    nz = int((scores != 0).sum())
    print(f"Saved {out}  shape={tuple(scores.shape)}  nonzero={nz}  "
          f"range=[{scores.min():.3f},{scores.max():.3f}]")


if __name__ == "__main__":
    main()
