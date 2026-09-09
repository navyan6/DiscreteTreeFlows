#!/usr/bin/env python3
"""
Table 6 (antibody affinity) metric helpers.

Paper metrics on Rodriguez Track C / OAS Ab trees:
  - Coverage@K at absolute Hamming radius e (default K=100, e∈{1,2,3,5})
  - SHM load error (|mean root→leaf edit frac gen − gt|)
  - CDR mutation recall (GT CDR substitutions recovered in any gen leaf)
  - Terminal diversity error (|mean pairwise leaf Hamming gen − gt|)
  - Lineage RF (NaN when topology forced or gen Newick missing)

CDR masks: prefer Rodriguez PCP codon coords; fallback IMGT-scaled stub.
"""

from __future__ import annotations

from typing import Iterable, Mapping, Sequence

__all__ = [
    "cdr_mask_imgt_stub",
    "cdr_mask_from_pcp_codon_coords",
    "cdr_mut_recall",
    "shm_load_error",
    "terminal_diversity_error",
    "coverage_at_k_abs",
    "mean_pairwise_hamming",
    "mean_shm_load",
    # backwards-compatible aliases used by eval_ab_maturation.py
    "cdr_mut_recall_stub",
    "shm_load_error_stub",
    "terminal_diversity_error_stub",
    "coverage_at_100_stub",
]


def cdr_mask_imgt_stub(seq_len: int, scheme: str = "imgt") -> list[bool]:
    """Approximate VH IMGT CDR ranges scaled to seq_len (AA). Fallback only."""
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


def cdr_mask_from_pcp_codon_coords(
    seq_len_aa: int,
    codon_coords: Mapping[str, tuple[int, int]] | None,
    *,
    nt_indexed: bool = True,
) -> list[bool]:
    """
    Build AA CDR mask from PCP ``cdr{1,2,3}_codon_{start,end}_heavy``.

    Rodriguez PCP stores NT indices of codon starts (inclusive). AA span for
    (start, end) is ``[start//3, end//3 + 1)`` when ``nt_indexed=True``.
    """
    if seq_len_aa <= 0:
        return []
    if not codon_coords:
        return cdr_mask_imgt_stub(seq_len_aa)
    mask = [False] * seq_len_aa
    for key in ("cdr1", "cdr2", "cdr3"):
        if key not in codon_coords:
            continue
        start, end = codon_coords[key]
        if start is None or end is None:
            continue
        try:
            start_i, end_i = int(start), int(end)
        except (TypeError, ValueError):
            continue
        if nt_indexed:
            lo = start_i // 3
            hi = end_i // 3 + 1
        else:
            lo, hi = start_i, end_i + 1
        lo = max(0, lo)
        hi = min(seq_len_aa, hi)
        for i in range(lo, hi):
            mask[i] = True
    if not any(mask):
        return cdr_mask_imgt_stub(seq_len_aa)
    return mask


def _mut_set(root: str, leaf: str, mask: Sequence[bool] | None = None) -> set[tuple[int, str]]:
    out: set[tuple[int, str]] = set()
    n = min(len(root), len(leaf))
    for i in range(n):
        if leaf[i] in "-." or root[i] in "-.":
            continue
        if mask is not None and (i >= len(mask) or not mask[i]):
            continue
        if leaf[i] != root[i]:
            out.add((i, leaf[i]))
    return out


def cdr_mut_recall(
    root: str,
    generated_leaves: Sequence[str],
    gt_leaves: Sequence[str],
    cdr_mask: Sequence[bool] | None = None,
) -> float:
    """Fraction of GT CDR substitutions (vs root) recovered in any generated leaf."""
    if not root or not gt_leaves:
        return float("nan")
    mask = list(cdr_mask) if cdr_mask is not None else cdr_mask_imgt_stub(len(root))
    if len(mask) != len(root):
        mask = cdr_mask_imgt_stub(len(root))

    gt_muts: set[tuple[int, str]] = set()
    for leaf in gt_leaves:
        gt_muts |= _mut_set(root, leaf, mask)
    if not gt_muts:
        return float("nan")

    recovered = 0
    for pos, aa in gt_muts:
        if any(pos < len(g) and g[pos] == aa for g in generated_leaves):
            recovered += 1
    return recovered / len(gt_muts)


def mean_shm_load(root: str, leaves: Sequence[str]) -> float:
    """Mean fractional Hamming distance root→leaf (gaps ignored in numerator)."""
    if not leaves or not root:
        return float("nan")
    loads = []
    for leaf in leaves:
        n = min(len(root), len(leaf))
        if n == 0:
            continue
        d = sum(
            1
            for i in range(n)
            if root[i] != leaf[i] and leaf[i] not in "-." and root[i] not in "-."
        )
        loads.append(d / n)
    return sum(loads) / len(loads) if loads else float("nan")


def shm_load_error(
    root: str,
    generated_leaves: Sequence[str],
    gt_leaves: Sequence[str],
) -> float:
    """|mean edit-fraction(root→gen) − mean edit-fraction(root→gt)|. Lower better."""
    a = mean_shm_load(root, generated_leaves)
    b = mean_shm_load(root, gt_leaves)
    if a != a or b != b:  # NaN
        return float("nan")
    return abs(a - b)


def mean_pairwise_hamming(leaves: Sequence[str]) -> float:
    if len(leaves) < 2:
        return float("nan")
    total = n = 0
    for i in range(len(leaves)):
        for j in range(i + 1, len(leaves)):
            a, b = leaves[i], leaves[j]
            m = min(len(a), len(b))
            if m == 0:
                continue
            total += sum(1 for k in range(m) if a[k] != b[k] and a[k] not in "-." and b[k] not in "-.")
            n += 1
    return total / n if n else float("nan")


def terminal_diversity_error(leaves_a: Sequence[str], leaves_b: Sequence[str]) -> float:
    """
    Absolute difference of mean pairwise (absolute) Hamming between leaf sets.

    Companion to antibody_benchmark leaf-diversity W1 (distributional); this is
    the scalar mean-error form for Table 6 'terminal diversity error'.
    """
    a = mean_pairwise_hamming(leaves_a)
    b = mean_pairwise_hamming(leaves_b)
    if a != a or b != b:
        return float("nan")
    return abs(a - b)


def coverage_at_k_abs(
    gt_leaves: Sequence[str],
    gen_leaves: Sequence[str],
    *,
    e: int = 2,
    k: int = 100,
    seed: int = 0,
) -> float:
    """
    Obs→gen Coverage@K at absolute Hamming radius e.

    Pool = first ``k`` gen leaves after deterministic shuffle (seed).
    """
    if not gt_leaves or not gen_leaves:
        return float("nan")
    pool = list(gen_leaves)
    if len(pool) > k:
        import random

        rng = random.Random(seed)
        rng.shuffle(pool)
        pool = pool[:k]
    try:
        from benchmarks.metrics.sequences import coverage_at_e

        return float(coverage_at_e(list(gt_leaves), pool, e=int(e)))
    except Exception:
        def ham(a: str, b: str) -> int:
            m = min(len(a), len(b))
            return sum(1 for i in range(m) if a[i] != b[i])

        hit = 0
        for g in gt_leaves:
            if any(ham(g, x) <= e for x in pool):
                hit += 1
        return hit / len(gt_leaves)


# --- aliases for older harness ---
cdr_mut_recall_stub = cdr_mut_recall
shm_load_error_stub = shm_load_error
terminal_diversity_error_stub = terminal_diversity_error


def coverage_at_100_stub(
    gt_leaves: Sequence[str],
    gen_leaves: Sequence[str],
    eps_frac: float = 0.02,
) -> float:
    """Legacy fractional Coverage@100; prefer coverage_at_k_abs for Table 6."""
    if not gt_leaves or not gen_leaves:
        return 0.0
    pool = list(gen_leaves)[:100]
    try:
        from benchmarks.metrics import sequences as S

        return float(S.coverage_at_k(list(gt_leaves), pool, eps_frac=eps_frac))
    except Exception:
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
    print("cdr_mut_recall", cdr_mut_recall(root, gen, gt, mask))
    print("shm_load_error", shm_load_error(root, gen, gt))
    print("terminal_diversity_error", terminal_diversity_error(gen, gt))
    print("coverage_at_k_abs e=2", coverage_at_k_abs(gt, gen + [gt[1]] * 98, e=2, k=100))
