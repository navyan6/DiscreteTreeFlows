"""
Sequence-recovery metrics for forecasting: compare a generated leaf-sequence set
to observed/true sequences, anchored on the root.

Used by both synthetic (Track A) and real viral blind forecasting (Track B1),
where generated leaves do NOT share identities with observed leaves — so metrics
are coverage / recovery based (best-of-K, coverage@K, mutation P/R/F1), not RF.

## Leaf vs tree (important)

Primary KPIs ``mut_recovery`` / ``site_recall`` / ``aa_acc_given_hit`` /
``cons_retention`` from ``positional_recovery`` are **leaf-only**:

  - GT: terminal (leaf) sequences vs root
  - Gen: best-matching generated **leaf** vs that GT leaf
  - Internal / ancestral nodes are **not** scored

This is intentional for forecasting (observed tips). Timing of when a mutation
appears on an internal edge does not enter the primary metric — only the leaf
AA matters. If GT mutations are tip-restricted while the model mutates on
internal edges (or vice versa), leaf scoring is still the fair end-state check;
use the tree-wide helpers below when you also want path-aggregated credit.

Tree-wide / fairer companions (do **not** replace the primary leaf metrics):
  - ``any_descendant_mut_recovery``: credit a GT leaf mut if **any** gen leaf
    reaches the GT AA (or mutates the site)
  - ``path_union_mutation_recovery``: union of mut sites over all GT leaves vs
    union over gen leaves (set recovery)
  - ``mutation_pr_f1`` (already leaf-union): precision/recall of mut sets
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter

import numpy as np

__all__ = [
    "hamming", "identity", "mutations_vs_root", "sites_vs_root",
    "best_of_k_identity", "min_hamming", "coverage_at_k",
    "coverage_at_e", "frac_gen_within_e",
    "mutation_pr_f1", "unique_mutations_recovered",
    "sitewise_entropy", "mutation_spectrum", "pairwise_distance_distribution",
    "align_index_map", "positional_recovery",
    "any_descendant_mut_recovery", "path_union_mutation_recovery",
]


# ── root-frame alignment ────────────────────────────────────────────────────
#
# Column-indexed metrics assume root[i], gt[i] and gen[i] are the same residue.
# That holds for generated leaves, which inherit the root's length because the
# models only substitute, but not for observed leaves, which carry their own
# indels. Where an observed leaf has a deletion, every residue downstream sits at
# a lower index than its homolog on the root and the comparison silently reads
# the wrong residue. See results/voc_threat_panel/SPIKE_FRAME_BUG.md.
#
# Exposure measured by scripts/audit_frame_all_datasets.py: 1.8% of HIV Env
# leaves share their root's frame, 47% for SARS-CoV-2 Spike, 74% for H1N1,
# 94% for H3N2.

_ALIGN_CACHE: dict[tuple[str, str], tuple[int, ...]] = {}
_ALIGN_CACHE_MAX = 20_000


def _digest(s: str) -> str:
    return hashlib.blake2b(s.encode(), digest_size=16).hexdigest()


def align_index_map(seq: str, root: str) -> tuple[int, ...]:
    """Map each index of ``root`` to the homologous index of ``seq``, or -1.

    Equal length means no indel relative to the root, so the identity map is
    exact and alignment is skipped — the common case, and what keeps this cheap
    on datasets like H3N2 where almost nothing has an indel.
    """
    if len(seq) == len(root):
        return tuple(range(len(root)))

    key = (_digest(root), _digest(seq))
    cached = _ALIGN_CACHE.get(key)
    if cached is not None:
        return cached

    mapping = [-1] * len(root)
    try:
        from Bio import Align

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -0.5
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        aln = aligner.align(root, seq)[0]
        for (rs, re_), (qs, _qe) in zip(*aln.aligned):
            for k in range(re_ - rs):
                mapping[rs + k] = qs + k
    except Exception:  # noqa: BLE001 - no Bio, or degenerate seq; truncate
        for i in range(min(len(root), len(seq))):
            mapping[i] = i

    out = tuple(mapping)
    if len(_ALIGN_CACHE) < _ALIGN_CACHE_MAX:
        _ALIGN_CACHE[key] = out
    return out


def _truncating_map(seq: str, root: str) -> tuple[int, ...]:
    """Legacy behaviour: index-for-index, truncated to the shorter sequence."""
    n = min(len(root), len(seq))
    return tuple(i if i < n else -1 for i in range(len(root)))


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


def coverage_at_e(targets: list[str], gen_seqs: list[str], e: int = 0) -> float:
    """
    Obs→gen coverage at absolute Hamming radius e:
      |{t ∈ targets : ∃g, d_H(t,g) ≤ e}| / |targets|
    """
    if not targets:
        return float("nan")
    if e < 0:
        raise ValueError("e must be >= 0")
    covered = sum(1 for t in targets if min_hamming(t, gen_seqs) <= e)
    return covered / len(targets)


def frac_gen_within_e(targets: list[str], gen_seqs: list[str], e: int = 0) -> float:
    """
    Gen→obs fraction at absolute Hamming radius e:
      |{g ∈ gen_seqs : ∃t, d_H(t,g) ≤ e}| / |gen_seqs|
    """
    if not gen_seqs:
        return float("nan")
    if e < 0:
        raise ValueError("e must be >= 0")
    if not targets:
        return 0.0
    near = sum(1 for g in gen_seqs if min_hamming(g, targets) <= e)
    return near / len(gen_seqs)


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


def positional_recovery(root: str, gt: str, gen: str, align: bool = True) -> dict:
    """
    **Primary leaf metric.** Split positions by root-vs-GT leaf and score one
    generated leaf:
      conserved (root==gt): model should keep root AA -> retention
      mutating  (root!=gt): model should reach GT AA  -> recovery

    Factorization of mut_recovery (exact when a wrong AA still counts as a site hit):
      mut_recovery      = P(gen==GT | root!=GT)
      site_recall       = P(gen!=root | root!=GT)          # mut-site hit rate
      aa_acc_given_hit  = P(gen==GT | root!=GT & gen!=root)
      => mut_recovery = site_recall * aa_acc_given_hit
         (when site_hits>0; if site_hits==0 then mut_recovery==0 and aa_acc is nan)

    Optional site_precision = P(root!=GT | gen!=root): among sites the model
    mutated away from root, fraction that were true mutating sites.

    Positions are scored in the **root's** coordinate frame. With ``align=True``
    (default) the GT and generated leaves are mapped onto that frame by pairwise
    alignment, so an indel-carrying leaf is read at its homologous residue rather
    than a shifted one; positions deleted on either leaf are dropped from both
    numerator and denominator and counted in ``n_skipped_indel``. Sequences of
    equal length skip alignment, so the fast path is unchanged and free.

    ``align=False`` restores the previous index-for-index behaviour, truncated to
    the shortest of the three. Kept only for reproducing pre-fix numbers.

    Does **not** look at internal nodes — see module docstring.
    """
    ix = align_index_map if align else _truncating_map
    gt_ix, gen_ix = ix(gt, root), ix(gen, root)

    mut_correct = mut_total = cons_correct = cons_total = 0
    site_hits = site_hit_correct = 0  # true mut sites where gen!=root; among those gen==GT
    gen_mut_total = gen_mut_true = 0  # gen!=root; among those root!=GT
    skipped = 0
    for i in range(len(root)):
        j, k = gt_ix[i], gen_ix[i]
        if j < 0 or k < 0:
            skipped += 1
            continue
        r, g, m = root[i], gt[j], gen[k]
        if r == g:
            cons_total += 1
            cons_correct += (m == r)
        else:
            mut_total += 1
            mut_correct += (m == g)
            if m != r:
                site_hits += 1
                site_hit_correct += (m == g)
        if m != r:
            gen_mut_total += 1
            gen_mut_true += (r != g)
    return {
        "mut_recovery": mut_correct / mut_total if mut_total else float("nan"),
        "cons_retention": cons_correct / cons_total if cons_total else float("nan"),
        "site_recall": site_hits / mut_total if mut_total else float("nan"),
        "aa_acc_given_hit": (
            site_hit_correct / site_hits if site_hits else float("nan")
        ),
        "site_precision": (
            gen_mut_true / gen_mut_total if gen_mut_total else float("nan")
        ),
        "mut_total": mut_total,
        "cons_total": cons_total,
        "site_hits": site_hits,
        "gen_mut_total": gen_mut_total,
        "n_scored": mut_total + cons_total,
        "n_skipped_indel": skipped,
    }


def any_descendant_mut_recovery(
    root: str,
    gt: str,
    gen_seqs: list[str],
    align: bool = True,
) -> dict:
    """
    Tree-wide companion to leaf ``positional_recovery`` for one GT leaf.

    For each site where root≠GT, credit recovery if **any** generated leaf has
    gen==GT at that site; site hit if any gen≠root. Conserved retention requires
    **all** gen leaves to keep root (strict) — also report soft mean retention.

    Scored in the root's frame with the same alignment handling as
    ``positional_recovery``; see that docstring for ``align``.

    Primary leaf ``mut_recovery`` is unchanged; this is an additional KPI.
    """
    if not gen_seqs:
        return {
            "mut_recovery_any_descendant": float("nan"),
            "site_recall_any_descendant": float("nan"),
            "cons_retention_all_gen": float("nan"),
            "mut_total": 0,
            "cons_total": 0,
            "n_skipped_indel": 0,
        }
    ix = align_index_map if align else _truncating_map
    gt_ix = ix(gt, root)
    gen_ix = [ix(s, root) for s in gen_seqs]

    mut_correct = mut_total = site_hits = 0
    cons_correct = cons_total = 0
    skipped = 0
    for i in range(len(root)):
        j = gt_ix[i]
        if j < 0 or any(m[i] < 0 for m in gen_ix):
            skipped += 1
            continue
        r, g = root[i], gt[j]
        gens = [s[m[i]] for s, m in zip(gen_seqs, gen_ix)]
        if r == g:
            cons_total += 1
            cons_correct += int(all(m == r for m in gens))
        else:
            mut_total += 1
            mut_correct += int(any(m == g for m in gens))
            site_hits += int(any(m != r for m in gens))
    return {
        "mut_recovery_any_descendant": (
            mut_correct / mut_total if mut_total else float("nan")
        ),
        "site_recall_any_descendant": (
            site_hits / mut_total if mut_total else float("nan")
        ),
        "cons_retention_all_gen": (
            cons_correct / cons_total if cons_total else float("nan")
        ),
        "mut_total": mut_total,
        "cons_total": cons_total,
        "n_skipped_indel": skipped,
    }


def path_union_mutation_recovery(
    root: str,
    gt_seqs: list[str],
    gen_seqs: list[str],
) -> dict:
    """
    Path-/tree-aggregated mutation set recovery (leaf unions vs root).

    GT mut sites = ∪_leaves sites_vs_root(root, gt_leaf)
    Gen mut sites = ∪_leaves sites_vs_root(root, gen_leaf)
    Gen substitutions = ∪_leaves mutations_vs_root(root, gen_leaf)

    Reports site-level recall/precision/F1 plus substitution-level recall of
    (pos, aa) pairs that appear on any GT leaf. Complements leaf-paired
    ``positional_recovery``; does not replace it.
    """
    site = mutation_pr_f1(gen_seqs, gt_seqs, root, level="site")
    sub = mutation_pr_f1(gen_seqs, gt_seqs, root, level="substitution")
    return {
        "mut_site_recall_path_union": site["recall"],
        "mut_site_precision_path_union": site["precision"],
        "mut_site_f1_path_union": site["f1"],
        "mut_sub_recall_path_union": sub["recall"],
        "mut_sub_precision_path_union": sub["precision"],
        "mut_sub_f1_path_union": sub["f1"],
        "n_gt_mut_sites": site["n_obs"],
        "n_gen_mut_sites": site["n_gen"],
        "n_gt_mut_subs": sub["n_obs"],
        "n_gen_mut_subs": sub["n_gen"],
    }
