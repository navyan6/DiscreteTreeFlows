"""
Distributional distances between a generated tree set (K samples) and a true
tree set (replicate simulations, or held-out reality).

Two families:
  * Summary-vector distances  — reduce each tree to a fixed feature vector, then
    compare the two clouds with per-feature Wasserstein-1, multivariate energy
    distance, and RBF-MMD.
  * Topology-distribution distances — tree-KL over rooted topologies and split-KL
    over clade (split) frequencies.

Pure numpy + scipy.stats.wasserstein_distance; operates on TreeState via
`benchmarks.metrics.trees`.
"""

from __future__ import annotations

import numpy as np

from src.tree_state import TreeState
from benchmarks.metrics import trees as T

__all__ = [
    "summarize_tree", "summary_matrix",
    "wasserstein_per_feature", "energy_distance", "mmd_rbf",
    "split_frequencies", "split_kl", "tree_kl",
    "distributional_report",
]

# Feature order for summary vectors (all cheap, model-free tree statistics).
SUMMARY_FEATURES = [
    "n_leaves", "n_nodes", "sackin", "colless", "cherries",
    "topo_height", "patristic_height", "width", "mean_branch_len",
]


def summarize_tree(tree: TreeState) -> np.ndarray:
    """Fixed-length feature vector for a tree (order = SUMMARY_FEATURES)."""
    bls = list(tree.branch_lengths.values())
    return np.array([
        float(tree.n_leaves()),
        float(tree.n_nodes()),
        float(T.sackin_index(tree)),
        float(T.colless_index(tree)),
        float(T.cherry_count(tree)),
        float(T.topological_height(tree)),
        float(T.patristic_height(tree)),
        float(T.tree_width(tree)),
        float(np.mean(bls)) if bls else 0.0,
    ], dtype=np.float64)


def summary_matrix(tree_list: list[TreeState]) -> np.ndarray:
    """[n_trees, n_features] stacked summary vectors."""
    if not tree_list:
        return np.empty((0, len(SUMMARY_FEATURES)))
    return np.vstack([summarize_tree(t) for t in tree_list])


# ── summary-vector distances ────────────────────────────────────────────────

def wasserstein_per_feature(X: np.ndarray, Y: np.ndarray) -> dict[str, float]:
    """1D Wasserstein-1 per feature between the two clouds."""
    from scipy.stats import wasserstein_distance
    out: dict[str, float] = {}
    for j, name in enumerate(SUMMARY_FEATURES):
        out[name] = float(wasserstein_distance(X[:, j], Y[:, j]))
    return out


def _standardize(X: np.ndarray, Y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Z-score both clouds by pooled mean/std so features are comparable."""
    pooled = np.vstack([X, Y])
    mu = pooled.mean(0)
    sd = pooled.std(0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd, (Y - mu) / sd


def _pairwise_euclidean(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    d2 = (
        (A * A).sum(1)[:, None]
        + (B * B).sum(1)[None, :]
        - 2.0 * A @ B.T
    )
    return np.sqrt(np.maximum(d2, 0.0))


def energy_distance(X: np.ndarray, Y: np.ndarray, standardize: bool = True) -> float:
    """
    Multivariate energy distance: 2 E|X-Y| - E|X-X'| - E|Y-Y'| (>= 0; 0 iff equal
    distributions in the limit).
    """
    if len(X) == 0 or len(Y) == 0:
        return float("nan")
    if standardize:
        X, Y = _standardize(X, Y)
    dxy = _pairwise_euclidean(X, Y).mean()
    dxx = _pairwise_euclidean(X, X).mean()
    dyy = _pairwise_euclidean(Y, Y).mean()
    return float(max(0.0, 2.0 * dxy - dxx - dyy))


def mmd_rbf(X: np.ndarray, Y: np.ndarray, gamma: float | None = None,
            standardize: bool = True) -> float:
    """
    Unbiased squared MMD with an RBF kernel. gamma defaults to the
    median-heuristic (1 / median pairwise squared distance).
    """
    if len(X) < 2 or len(Y) < 2:
        return float("nan")
    if standardize:
        X, Y = _standardize(X, Y)

    def sqdist(A, B):
        return np.maximum(
            (A * A).sum(1)[:, None] + (B * B).sum(1)[None, :] - 2.0 * A @ B.T, 0.0
        )

    dxx, dyy, dxy = sqdist(X, X), sqdist(Y, Y), sqdist(X, Y)
    if gamma is None:
        med = np.median(np.concatenate([dxx.ravel(), dyy.ravel(), dxy.ravel()]))
        gamma = 1.0 / med if med > 0 else 1.0

    kxx, kyy, kxy = np.exp(-gamma * dxx), np.exp(-gamma * dyy), np.exp(-gamma * dxy)
    m, n = len(X), len(Y)
    # unbiased: drop diagonal self-terms
    np.fill_diagonal(kxx, 0.0)
    np.fill_diagonal(kyy, 0.0)
    term_xx = kxx.sum() / (m * (m - 1))
    term_yy = kyy.sum() / (n * (n - 1))
    term_xy = kxy.mean()
    return float(max(0.0, term_xx + term_yy - 2.0 * term_xy))


# ── topology-distribution distances ─────────────────────────────────────────

def _topology_key(tree: TreeState) -> frozenset:
    """Canonical rooted-topology identity: the set of informative clades."""
    return frozenset(T.clades(tree))


def split_frequencies(tree_list: list[TreeState]) -> dict[frozenset, float]:
    """Fraction of trees containing each clade (split)."""
    n = len(tree_list)
    if n == 0:
        return {}
    counts: dict[frozenset, int] = {}
    for t in tree_list:
        for c in T.clades(t):
            counts[c] = counts.get(c, 0) + 1
    return {c: k / n for c, k in counts.items()}


def split_kl(gen: list[TreeState], true: list[TreeState], eps: float = 1e-3) -> float:
    """
    Sum over the union of clades of the Bernoulli KL between generated and true
    split-presence frequencies: sum_c KL( Bern(p_c) || Bern(q_c) ).
    """
    pf, qf = split_frequencies(gen), split_frequencies(true)
    clades = set(pf) | set(qf)
    total = 0.0
    for c in clades:
        p = min(1 - eps, max(eps, pf.get(c, 0.0)))
        q = min(1 - eps, max(eps, qf.get(c, 0.0)))
        total += p * np.log(p / q) + (1 - p) * np.log((1 - p) / (1 - q))
    return float(total)


def tree_kl(gen: list[TreeState], true: list[TreeState], eps: float = 1e-6) -> float:
    """
    KL( P_gen || P_true ) over the categorical distribution of rooted topologies
    (Laplace-smoothed over the union of observed topologies).
    """
    from collections import Counter
    pg, pt = Counter(_topology_key(t) for t in gen), Counter(_topology_key(t) for t in true)
    keys = set(pg) | set(pt)
    ng, nt = len(gen), len(true)
    V = len(keys)
    total = 0.0
    for k in keys:
        p = (pg.get(k, 0) + eps) / (ng + eps * V)
        q = (pt.get(k, 0) + eps) / (nt + eps * V)
        total += p * np.log(p / q)
    return float(total)


# ── convenience report ──────────────────────────────────────────────────────

def distributional_report(gen: list[TreeState], true: list[TreeState]) -> dict:
    """All distributional distances between a generated and a true tree set."""
    Xg, Xt = summary_matrix(gen), summary_matrix(true)
    return {
        "wasserstein": wasserstein_per_feature(Xg, Xt),
        "energy_distance": energy_distance(Xg, Xt),
        "mmd_rbf": mmd_rbf(Xg, Xt),
        "split_kl": split_kl(gen, true),
        "tree_kl": tree_kl(gen, true),
        "n_gen": len(gen),
        "n_true": len(true),
    }
