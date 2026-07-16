"""
Leaf matching + sequence-topology coupling.

For blind forecasting, generated leaves have no shared identities with observed
leaves, so we match them by sequence (Hungarian / min-cost assignment) and then
ask whether TreeSBM placed those sequences on the tree consistently with reality
— i.e. it did NOT generate realistic sequences and a realistic tree
*independently*. That coupling is measured via patristic-distance agreement.
"""

from __future__ import annotations

import numpy as np

from src.tree_state import TreeState
from benchmarks.metrics import trees as T
from benchmarks.metrics.sequences import hamming

__all__ = [
    "patristic_matrix", "seq_patristic_correlation",
    "hungarian_match", "matched_patristic_agreement",
]


def _lca_getter(tree: TreeState):
    pm = T.parent_map(tree)
    depth = T.node_depths(tree)

    def lca(u: str, v: str) -> str:
        while depth[u] > depth[v]:
            u = pm[u]
        while depth[v] > depth[u]:
            v = pm[v]
        while u != v:
            u, v = pm[u], pm[v]
        return u
    return lca


def patristic_matrix(tree: TreeState, leaves: list[str] | None = None
                     ) -> tuple[list[str], np.ndarray]:
    """Pairwise patristic (branch-length) distances among the given leaves."""
    if leaves is None:
        leaves = T.leaf_labels(tree)
    times = T.node_times(tree)
    lca = _lca_getter(tree)
    n = len(leaves)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a = lca(leaves[i], leaves[j])
            d = times[leaves[i]] + times[leaves[j]] - 2.0 * times[a]
            M[i, j] = M[j, i] = d
    return leaves, M


def _upper(M: np.ndarray) -> np.ndarray:
    iu = np.triu_indices(M.shape[0], k=1)
    return M[iu]


def seq_patristic_correlation(tree: TreeState,
                              node_seqs: dict[str, str] | None = None,
                              max_leaves: int = 120, seed: int = 0) -> float:
    """
    Pearson correlation between pairwise sequence distance and pairwise patristic
    distance among a tree's leaves (within-tree sequence-topology coupling).
    """
    import random
    seqs = node_seqs if node_seqs is not None else tree.node_seqs
    leaves = T.leaf_labels(tree)
    if len(leaves) > max_leaves:
        leaves = random.Random(seed).sample(leaves, max_leaves)
    if len(leaves) < 3:
        return float("nan")
    _, P = patristic_matrix(tree, leaves)
    n = len(leaves)
    Sq = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            Sq[i, j] = Sq[j, i] = hamming(seqs[leaves[i]], seqs[leaves[j]])
    pv, sv = _upper(P), _upper(Sq)
    if pv.std() == 0 or sv.std() == 0:
        return float("nan")
    return float(np.corrcoef(pv, sv)[0, 1])


def hungarian_match(true_seqs: list[str], gen_seqs: list[str]) -> dict:
    """
    Min-Hamming assignment of true sequences to generated sequences.
    Handles rectangular inputs (matches min(n_true, n_gen)).
    """
    from scipy.optimize import linear_sum_assignment
    n, m = len(true_seqs), len(gen_seqs)
    cost = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            cost[i, j] = hamming(true_seqs[i], gen_seqs[j])
    row, col = linear_sum_assignment(cost)
    matched = [(int(r), int(c)) for r, c in zip(row, col)]
    ham = [cost[r, c] for r, c in matched]
    lens = [min(len(true_seqs[r]), len(gen_seqs[c])) for r, c in matched]
    ident = [1 - h / L if L else 0.0 for h, L in zip(ham, lens)]
    return {
        "matches": matched,
        "mean_matched_hamming": float(np.mean(ham)) if ham else float("nan"),
        "mean_matched_identity": float(np.mean(ident)) if ident else float("nan"),
        "n_matched": len(matched),
    }


def matched_patristic_agreement(true_tree: TreeState, gen_tree: TreeState,
                                true_seqs: dict[str, str] | None = None,
                                gen_seqs: dict[str, str] | None = None,
                                max_leaves: int = 100, seed: int = 0) -> dict:
    """
    Match true leaves to generated leaves by sequence, then compare the two
    patristic-distance matrices over the matched leaves.

    Returns Pearson correlation and normalized RMSE — high correlation means the
    generated sequence placement is consistent with the true tree geometry
    (sequence and topology were generated jointly, not independently).
    """
    import random
    ts = true_seqs if true_seqs is not None else true_tree.node_seqs
    gs = gen_seqs if gen_seqs is not None else gen_tree.node_seqs

    t_leaves = T.leaf_labels(true_tree)
    g_leaves = T.leaf_labels(gen_tree)
    if len(t_leaves) > max_leaves:
        t_leaves = random.Random(seed).sample(t_leaves, max_leaves)

    match = hungarian_match([ts[l] for l in t_leaves], [gs[l] for l in g_leaves])
    pairs = match["matches"]
    if len(pairs) < 3:
        return {"patristic_corr": float("nan"), "patristic_nrmse": float("nan"),
                **match}

    t_matched = [t_leaves[r] for r, _ in pairs]
    g_matched = [g_leaves[c] for _, c in pairs]
    _, Pt = patristic_matrix(true_tree, t_matched)
    _, Pg = patristic_matrix(gen_tree, g_matched)
    pt, pg = _upper(Pt), _upper(Pg)

    corr = (float(np.corrcoef(pt, pg)[0, 1])
            if pt.std() > 0 and pg.std() > 0 else float("nan"))
    scale = pt.max() if pt.max() > 0 else 1.0
    nrmse = float(np.sqrt(np.mean((pt - pg) ** 2)) / scale)
    return {"patristic_corr": corr, "patristic_nrmse": nrmse, **match}
