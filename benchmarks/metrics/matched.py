"""
Sequence-matched topology metrics.

Generated leaves have no shared identities with the reference leaves, so before
any RF/quartet computation we match generated terminals to reference terminals
by minimum-cost bipartite assignment on normalized amino-acid Hamming distance,
relabel the generated leaves accordingly, then compute normalized RF / quartet.
Also provides the matched terminal-sequence edit distance.
"""

from __future__ import annotations

import numpy as np

from src.tree_state import TreeState
from benchmarks.metrics import trees as T
from benchmarks.metrics.sequences import hamming, identity

__all__ = [
    "match_leaves", "relabel_leaves", "sequence_matched_rf",
    "quartet_distance", "terminal_edit_distance",
]


def _norm_hamming(a: str, b: str) -> float:
    L = min(len(a), len(b))
    return hamming(a, b) / L if L else 1.0


def match_leaves(gen: TreeState, ref: TreeState) -> dict[str, str]:
    """
    Min-cost matching of generated leaves -> reference leaves by normalized AA
    Hamming. Requires equal leaf counts. Returns {gen_leaf: ref_leaf}.
    """
    from scipy.optimize import linear_sum_assignment
    g_leaves, r_leaves = T.leaf_labels(gen), T.leaf_labels(ref)
    if len(g_leaves) != len(r_leaves):
        raise ValueError(f"leaf-count mismatch: gen {len(g_leaves)} vs ref {len(r_leaves)}")
    cost = np.array([[_norm_hamming(gen.node_seqs[g], ref.node_seqs[r])
                      for r in r_leaves] for g in g_leaves])
    row, col = linear_sum_assignment(cost)
    return {g_leaves[i]: r_leaves[j] for i, j in zip(row, col)}


def relabel_leaves(tree: TreeState, mapping: dict[str, str]) -> TreeState:
    """Return a copy of `tree` with leaves renamed per `mapping` (internal nodes unchanged)."""
    def rn(n):
        return mapping.get(n, n)
    return TreeState(
        node_ids=[rn(n) for n in tree.node_ids],
        root_id=rn(tree.root_id),
        edges=[(rn(p), rn(c)) for p, c in tree.edges],
        branch_lengths={(rn(p), rn(c)): v for (p, c), v in tree.branch_lengths.items()},
        node_seqs={rn(n): s for n, s in tree.node_seqs.items()},
        active_leaves=[rn(n) for n in tree.active_leaves],
    )


def _matched(gen: TreeState, ref: TreeState) -> TreeState:
    return relabel_leaves(gen, match_leaves(gen, ref))


def sequence_matched_rf(gen: TreeState, ref: TreeState) -> float:
    """Normalized rooted RF after sequence-matched relabeling of generated leaves."""
    return T.normalized_rf(_matched(gen, ref), ref)


def quartet_distance(gen: TreeState, ref: TreeState) -> float:
    """
    Normalized quartet distance (tqDist) after sequence-matched relabeling.
    Requires the `tqdist` package (pip install tqdist). Raises if unavailable so
    the caller can record it as a missing dependency rather than a fake number.
    """
    try:
        from tqdist import quartet_distance as _qd
    except ImportError as e:  # pragma: no cover
        raise ImportError("quartet_distance needs `tqdist` (pip install tqdist)") from e
    g = _matched(gen, ref)
    # tqdist compares unrooted quartets over the shared leaf label set. Internal
    # node names must be OMITTED: different tree sources name internal nodes
    # differently (e.g. dendropy sims vs Bio.Phylo-parsed real trees), and tqdist
    # treats any embedded internal label as a taxon that must match between the
    # two trees, aborting on mismatch. Only leaves (already sequence-matched)
    # should be labeled.
    try:
        return float(_qd(T.to_newick(g, with_names=False), T.to_newick(ref, with_names=False)))
    except Exception:
        return float("nan")


def terminal_edit_distance(gen: TreeState, ref: TreeState) -> dict:
    """
    Mean matched normalized terminal-sequence distance (1 - matched identity),
    plus the matched pair distances for saving distributions.
    """
    mapping = match_leaves(gen, ref)
    dists = [_norm_hamming(gen.node_seqs[g], ref.node_seqs[r])
             for g, r in mapping.items()]
    return {
        "mean": float(np.mean(dists)) if dists else float("nan"),
        "matched_distances": dists,
    }
