"""
Sequence-recovery metrics for forecasting: compare a generated leaf-sequence set
to observed/true sequences, anchored on the root.

Used by both synthetic (Track A) and real viral blind forecasting (Track B1),
where generated leaves do NOT share identities with observed leaves — so metrics
are coverage / recovery based (best-of-K, coverage@K, mutation P/R/F1), not RF.
"""

from __future__ import annotations

import math
from collections import Counter

import numpy as np

__all__ = [
    "hamming", "identity", "mutations_vs_root", "sites_vs_root",
    "best_of_k_identity", "min_hamming", "coverage_at_k",
    "mutation_pr_f1", "unique_mutations_recovered",
    "sitewise_entropy", "mutation_spectrum", "pairwise_distance_distribution",
    "positional_recovery",
]


def hamming(a: str, b: str) -> int:
    L = min(len(a), len(b))
    return sum(1 for i in range(L) if a[i] != b[i])


def identity(a: str, b: str) -> float:
    L = min(len(a), len(b))
    return sum(a[i] == b[i] for i in range(L)) / L if L else 0.0


def mutations_vs_root(root: str, seq: str) -> set[tuple[int, str]]:
    """Set of (position, mutant_aa) where seq differs from root."""
    L = min(len(root), len(seq))
    return {(i, seq[i]) for i in range(L) if seq[i] != root[i]}


def sites_vs_root(root: str, seq: str) -> set[int]:
    L = min(len(root), len(seq))
    return {i for i in range(L) if seq[i] != root[i]}


# ── best-of-K / coverage ────────────────────────────────────────────────────

def best_of_k_identity(target: str, gen_seqs: list[str]) -> float:
    """Max sequence identity of any generated leaf to the target."""
    return max((identity(target, g) for g in gen_seqs), default=0.0)


def min_hamming(target: str, gen_seqs: list[str]) -> int:
    return min((hamming(target, g) for g in gen_seqs), default=len(target))


def coverage_at_k(targets: list[str], gen_seqs: list[str], eps_frac: float = 0.0) -> float:
    """
    Fraction of target sequences within `eps_frac` fractional-Hamming of at least
    one generated leaf. eps_frac=0 requires an exact-length match.
    """
    if not targets:
        return float("nan")
    covered = 0
    for t in targets:
        L = len(t)
        thresh = eps_frac * L
        if any(hamming(t, g) <= thresh for g in gen_seqs):
            covered += 1
    return covered / len(targets)


# ── mutation precision / recall / F1 ────────────────────────────────────────

def _mut_set(seqs: list[str], root: str, level: str) -> set:
    acc: set = set()
    for s in seqs:
        if level == "substitution":
            acc |= mutations_vs_root(root, s)
        elif level == "site":
            acc |= sites_vs_root(root, s)
        else:
            raise ValueError("level must be 'substitution' or 'site'")
    return acc


def mutation_pr_f1(gen_seqs: list[str], true_seqs: list[str], root: str,
                   level: str = "substitution") -> dict[str, float]:
    """
    Precision/recall/F1 of the generated mutation set against the observed
    (true future) mutation set, both taken as the union over leaves vs root.
    `level`: 'substitution' = (pos, aa) pairs; 'site' = positions only.
    """
    gen = _mut_set(gen_seqs, root, level)
    obs = _mut_set(true_seqs, root, level)
    tp = len(gen & obs)
    precision = tp / len(gen) if gen else 0.0
    recall = tp / len(obs) if obs else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {"precision": precision, "recall": recall, "f1": f1,
            "n_gen": len(gen), "n_obs": len(obs), "n_recovered": tp}


def unique_mutations_recovered(gen_seqs: list[str], true_seqs: list[str], root: str,
                               level: str = "substitution") -> int:
    """Count of distinct observed future mutations present in the generated set."""
    return len(_mut_set(gen_seqs, root, level) & _mut_set(true_seqs, root, level))


# ── distributional sequence summaries ───────────────────────────────────────

def sitewise_entropy(seqs: list[str]) -> np.ndarray:
    """Per-position Shannon entropy (nats) over a set of sequences."""
    if not seqs:
        return np.zeros(0)
    L = min(len(s) for s in seqs)
    ent = np.zeros(L)
    n = len(seqs)
    for i in range(L):
        counts = Counter(s[i] for s in seqs)
        ent[i] = -sum((c / n) * math.log(c / n) for c in counts.values())
    return ent


def mutation_spectrum(seqs: list[str], root: str) -> dict[str, int]:
    """Aggregate count of each substitution type 'X>Y' across leaves (vs root)."""
    spec: Counter = Counter()
    for s in seqs:
        L = min(len(root), len(s))
        for i in range(L):
            if s[i] != root[i]:
                spec[f"{root[i]}>{s[i]}"] += 1
    return dict(spec)


def pairwise_distance_distribution(seqs: list[str], max_pairs: int = 2000,
                                   seed: int = 0) -> list[float]:
    """Fractional-Hamming distances over a (subsampled) set of leaf pairs."""
    import random
    rng = random.Random(seed)
    pairs = [(i, j) for i in range(len(seqs)) for j in range(i + 1, len(seqs))]
    if len(pairs) > max_pairs:
        pairs = rng.sample(pairs, max_pairs)
    out = []
    for i, j in pairs:
        L = min(len(seqs[i]), len(seqs[j]))
        out.append(hamming(seqs[i], seqs[j]) / L if L else 0.0)
    return out


def positional_recovery(root: str, gt: str, gen: str) -> dict:
    """
    Split positions by root-vs-GT and score the generated sequence:
      conserved (root==gt): model should keep root AA -> retention
      mutating  (root!=gt): model should reach GT AA  -> recovery
    """
    L = min(len(root), len(gt), len(gen))
    mut_correct = mut_total = cons_correct = cons_total = 0
    for i in range(L):
        r, g, m = root[i], gt[i], gen[i]
        if r == g:
            cons_total += 1
            cons_correct += (m == r)
        else:
            mut_total += 1
            mut_correct += (m == g)
    return {
        "mut_recovery": mut_correct / mut_total if mut_total else float("nan"),
        "cons_retention": cons_correct / cons_total if cons_total else float("nan"),
        "mut_total": mut_total,
        "cons_total": cons_total,
    }
