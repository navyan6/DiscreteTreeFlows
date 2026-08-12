#!/usr/bin/env python3
"""
Table 5 / Table 6 (antibody) metric helper stubs.

Paper draft numbering: antibody affinity block is Table 6 in ICLR_TABLE_FILL_PLAN;
ops shorthand "Table 5 Ab baselines" uses the same metrics.

CDR mut. recall / SHM load error / terminal diversity — scaffolding only.
Wire these into scripts/eval_ab_maturation.py once clone trees + ANARCI CDR
masks exist. Do not invent paper numbers from these stubs.
"""

from __future__ import annotations

from typing import Sequence


def cdr_mask_imgt_stub(seq_len: int, scheme: str = "imgt") -> list[bool]:
    """
    Placeholder CDR mask (VH IMGT-ish ranges scaled to seq_len).

    Real path: ANARCI / IgBLAST CDR1–3 coordinates per sequence.
    Returns a boolean list length seq_len (True = CDR position).
    """
    # Classic IMGT VH CDR approx on 128-aligned: CDR1 27–38, CDR2 56–65, CDR3 105–117
    # Scale linearly for variable unaligned VH length (very rough).
    _ = scheme
    if seq_len <= 0:
        return []
    scale = seq_len / 128.0
    ranges = [(27, 38), (56, 65), (105, 117)]
    mask = [False] * seq_len
    for a, b in ranges:
        lo = max(0, int(a * scale) - 1)
        hi = min(seq_len, int(b * scale))
        for i in range(lo, hi):
            mask[i] = True
    return mask


def cdr_mut_recall_stub(
    root: str,
    generated_leaves: Sequence[str],
    gt_leaves: Sequence[str],
    cdr_mask: Sequence[bool] | None = None,
) -> float:
    """
    Fraction of GT CDR substitutions (vs root) recovered in any generated leaf.

    Stub: identity-position amino-acid compare; no alignment. Returns 0.0 if
    empty inputs. NOT a FINAL Table 6 number.
    """
    if not root or not gt_leaves:
        return 0.0
    L = len(root)
    mask = list(cdr_mask) if cdr_mask is not None else cdr_mask_imgt_stub(L)
    if len(mask) != L:
        mask = cdr_mask_imgt_stub(L)

    gt_muts: set[tuple[int, str]] = set()
    for leaf in gt_leaves:
        for i in range(min(L, len(leaf))):
            if mask[i] and leaf[i] != root[i] and leaf[i] not in "-.":
                gt_muts.add((i, leaf[i]))
    if not gt_muts:
        return 0.0

    recovered = 0
    for pos, aa in gt_muts:
        for gen in generated_leaves:
            if pos < len(gen) and gen[pos] == aa:
                recovered += 1
                break
    return recovered / len(gt_muts)


def shm_load_error_stub(
    root: str,
    generated_leaves: Sequence[str],
    gt_leaves: Sequence[str],
) -> float:
    """
    |mean edit-fraction(root→gen) − mean edit-fraction(root→gt)|.

    Stub Hamming on min length; gaps ignored. Lower is better. NOT FINAL.
    """
    def mean_load(leaves: Sequence[str]) -> float:
        if not leaves or not root:
            return 0.0
        loads = []
        for leaf in leaves:
            n = min(len(root), len(leaf))
            if n == 0:
                continue
            d = sum(1 for i in range(n) if root[i] != leaf[i] and leaf[i] not in "-.")
            loads.append(d / n)
        return sum(loads) / len(loads) if loads else 0.0

    return abs(mean_load(generated_leaves) - mean_load(gt_leaves))


def terminal_diversity_error_stub(leaves_a: Sequence[str], leaves_b: Sequence[str]) -> float:
    """
    Absolute difference of mean pairwise Hamming diversity between two leaf sets.

    Stub only — use for scaffolding eval_ab_maturation.py later.
    """
    def mean_pairwise(leaves: Sequence[str]) -> float:
        if len(leaves) < 2:
            return 0.0
        total = n = 0
        for i in range(len(leaves)):
            for j in range(i + 1, len(leaves)):
                a, b = leaves[i], leaves[j]
                m = min(len(a), len(b))
                if m == 0:
                    continue
                total += sum(1 for k in range(m) if a[k] != b[k]) / m
                n += 1
        return total / n if n else 0.0

    return abs(mean_pairwise(leaves_a) - mean_pairwise(leaves_b))


def coverage_at_100_stub(
    gt_leaves: Sequence[str],
    gen_leaves: Sequence[str],
    eps_frac: float = 0.02,
) -> float:
    """
    Coverage@100 proxy: fraction of GT leaves with a gen neighbor within eps.

    Prefers benchmarks.metrics.sequences.coverage_at_k when importable.
    """
    if not gt_leaves or not gen_leaves:
        return 0.0
    pool = list(gen_leaves)[:100]
    try:
        from benchmarks.metrics import sequences as S

        return float(S.coverage_at_k(list(gt_leaves), pool, eps_frac=eps_frac))
    except Exception:
        # Hamming identity fallback
        def ident(a: str, b: str) -> float:
            n = min(len(a), len(b))
            if n == 0:
                return 0.0
            return sum(1 for i in range(n) if a[i] == b[i]) / n

        hit = 0
        for g in gt_leaves:
            if any(ident(g, x) >= 1.0 - eps_frac for x in pool):
                hit += 1
        return hit / len(gt_leaves)


if __name__ == "__main__":
    root = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQG"
    gt = [root[:20] + "A" + root[21:], root[:30] + "C" + root[31:]]
    gen = [gt[0], root]
    mask = cdr_mask_imgt_stub(len(root))
    print("cdr_mut_recall_stub", cdr_mut_recall_stub(root, gen, gt, mask))
    print("shm_load_error_stub", shm_load_error_stub(root, gen, gt))
    print("terminal_diversity_error_stub", terminal_diversity_error_stub(gen, gt))
    print("coverage_at_100_stub", coverage_at_100_stub(gt, gen + [gt[1]] * 98))
