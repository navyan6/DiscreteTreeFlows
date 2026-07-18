"""
Branch-length Wasserstein distances between a generated tree and a reference.

Reports the primary aggregate (all branches) plus the components (internal,
pendant, root-to-tip) which the driver saves separately.
"""

from __future__ import annotations

from src.tree_state import TreeState
from benchmarks.metrics import trees as T

__all__ = ["branch_length_sets", "branch_length_wasserstein"]


def branch_length_sets(tree: TreeState) -> dict[str, list[float]]:
    leaves = set(T.leaf_labels(tree))
    all_bl = [float(v) for v in tree.branch_lengths.values()]
    pendant = [float(v) for (p, c), v in tree.branch_lengths.items() if c in leaves]
    internal = [float(v) for (p, c), v in tree.branch_lengths.items() if c not in leaves]
    times = T.node_times(tree)
    roottotip = [times[l] for l in leaves]
    return {"all": all_bl, "internal": internal, "pendant": pendant,
            "roottotip": roottotip}


def branch_length_wasserstein(gen: TreeState, ref: TreeState) -> dict[str, float]:
    """W1 per component; primary table value is the 'all' entry."""
    from scipy.stats import wasserstein_distance
    g, r = branch_length_sets(gen), branch_length_sets(ref)
    out = {}
    for k in ("all", "internal", "pendant", "roottotip"):
        gv, rv = g[k], r[k]
        out[k] = float(wasserstein_distance(gv, rv)) if gv and rv else float("nan")
    return out
