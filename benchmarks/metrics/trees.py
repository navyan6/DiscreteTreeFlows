"""
Tree topology / shape metrics and newick I/O for TreeState.

All metrics are pure-Python and operate directly on `src.tree_state.TreeState`
(rooted, bifurcating trees with parent->child edges and branch lengths).

Topology comparison uses rooted **clades** (leaf-descendant sets of internal
nodes), which is the natural object for the rooted trees TreeSBM generates and
also yields split frequencies for split-KL in `distributions.py`.

`ete3` is optional (only `to_ete3`); everything else is dependency-light
(stdlib + the TreeState API), so this module imports and tests without ete3.
"""

from __future__ import annotations

import re
from collections import defaultdict, deque

from src.tree_state import TreeState

__all__ = [
    "children_map", "parent_map", "leaf_labels", "node_depths", "node_times",
    "descendant_leaves", "clades", "rf_distance", "normalized_rf",
    "sackin_index", "colless_index", "cherry_count",
    "topological_height", "patristic_height", "tree_width", "ltt",
    "to_newick", "to_ete3",
]


# ── topology helpers ────────────────────────────────────────────────────────

def children_map(tree: TreeState) -> dict[str, list[str]]:
    cm: dict[str, list[str]] = defaultdict(list)
    for p, c in tree.edges:
        cm[p].append(c)
    return dict(cm)


def parent_map(tree: TreeState) -> dict[str, str]:
    return {c: p for p, c in tree.edges}


def leaf_labels(tree: TreeState) -> list[str]:
    cm = children_map(tree)
    return [n for n in tree.node_ids if n not in cm]


def node_depths(tree: TreeState) -> dict[str, int]:
    """Topological depth (edges from root) per node, via BFS."""
    cm = children_map(tree)
    depths = {tree.root_id: 0}
    q = deque([tree.root_id])
    while q:
        n = q.popleft()
        for c in cm.get(n, []):
            depths[c] = depths[n] + 1
            q.append(c)
    return depths


def node_times(tree: TreeState) -> dict[str, float]:
    """Time (cumulative branch length from root) per node, via BFS."""
    cm = children_map(tree)
    times = {tree.root_id: 0.0}
    q = deque([tree.root_id])
    while q:
        n = q.popleft()
        for c in cm.get(n, []):
            times[c] = times[n] + float(tree.branch_lengths.get((n, c), 0.0))
            q.append(c)
    return times


def descendant_leaves(tree: TreeState) -> dict[str, frozenset[str]]:
    """Leaf-descendant set for every node (post-order accumulation)."""
    cm = children_map(tree)
    leaves = set(leaf_labels(tree))
    out: dict[str, frozenset[str]] = {}

    # iterative post-order to avoid recursion limits on deep/caterpillar trees
    order: list[str] = []
    stack = [tree.root_id]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(cm.get(n, []))
    for n in reversed(order):
        if n in leaves:
            out[n] = frozenset([n])
        else:
            acc: set[str] = set()
            for c in cm.get(n, []):
                acc |= out[c]
            out[n] = frozenset(acc)
    return out


# ── topology comparison (clade / RF) ────────────────────────────────────────

def clades(tree: TreeState, min_size: int = 2) -> set[frozenset[str]]:
    """
    Rooted clades: leaf-descendant sets of size >= min_size, excluding the
    full leaf set (uninformative). These are the rooted analogue of splits.
    """
    dl = descendant_leaves(tree)
    all_leaves = dl[tree.root_id]
    return {
        s for s in dl.values()
        if min_size <= len(s) < len(all_leaves)
    }


def rf_distance(t1: TreeState, t2: TreeState) -> int:
    """
    Rooted Robinson-Foulds: size of the symmetric difference of clade sets.
    Requires identical leaf label sets (raises otherwise).
    """
    l1, l2 = set(leaf_labels(t1)), set(leaf_labels(t2))
    if l1 != l2:
        raise ValueError(
            f"RF needs matching leaf sets; got |t1|={len(l1)}, |t2|={len(l2)}, "
            f"shared={len(l1 & l2)}"
        )
    c1, c2 = clades(t1), clades(t2)
    return len(c1 ^ c2)


def normalized_rf(t1: TreeState, t2: TreeState) -> float:
    """RF normalized by the maximum possible (|c1| + |c2|); 0 = identical."""
    c1, c2 = clades(t1), clades(t2)
    denom = len(c1) + len(c2)
    return len(c1 ^ c2) / denom if denom else 0.0


# ── shape / balance statistics ──────────────────────────────────────────────

def sackin_index(tree: TreeState) -> int:
    """Sum of leaf topological depths (higher = more imbalanced)."""
    d = node_depths(tree)
    return sum(d[l] for l in leaf_labels(tree))


def colless_index(tree: TreeState) -> int:
    """
    Sum over internal nodes of |L - R| leaf counts (bifurcating trees).
    Multifurcations contribute max-min across children.
    """
    cm = children_map(tree)
    dl = descendant_leaves(tree)
    total = 0
    for n, kids in cm.items():
        sizes = [len(dl[c]) for c in kids]
        if len(sizes) >= 2:
            total += max(sizes) - min(sizes)
    return total


def cherry_count(tree: TreeState) -> int:
    """Number of internal nodes whose (2) children are both leaves."""
    cm = children_map(tree)
    leaves = set(leaf_labels(tree))
    return sum(1 for kids in cm.values()
               if len(kids) == 2 and all(c in leaves for c in kids))


def topological_height(tree: TreeState) -> int:
    d = node_depths(tree)
    return max(d.values()) if d else 0


def patristic_height(tree: TreeState) -> float:
    """Max root-to-leaf branch-length distance."""
    t = node_times(tree)
    leaves = leaf_labels(tree)
    return max((t[l] for l in leaves), default=0.0)


def tree_width(tree: TreeState) -> int:
    """Max number of nodes at any single topological depth."""
    d = node_depths(tree)
    per_level: dict[int, int] = defaultdict(int)
    for depth in d.values():
        per_level[depth] += 1
    return max(per_level.values()) if per_level else 0


def ltt(tree: TreeState, n_points: int = 100) -> tuple[list[float], list[int]]:
    """
    Lineages-through-time on a uniform time grid over [0, patristic_height].

    Lineages at time t = number of edges (p, c) whose time interval
    [time(p), time(c)] contains t (standard reconstructed-tree LTT that
    handles serially-sampled tips).
    """
    times = node_times(tree)
    intervals = [(times[p], times[c]) for p, c in tree.edges]
    h = patristic_height(tree)
    if h <= 0 or not intervals:
        return [0.0], [max(1, tree.n_leaves())]
    grid = [h * i / (n_points - 1) for i in range(n_points)]
    counts = [sum(1 for a, b in intervals if a <= t < b or (t == h and b == h))
              for t in grid]
    return grid, counts


# ── newick / ete3 I/O ───────────────────────────────────────────────────────

_UNSAFE = re.compile(r"[(),:;\[\]\s]")


def _safe(name: str) -> str:
    return _UNSAFE.sub("_", name)


def to_newick(tree: TreeState, with_names: bool = True) -> str:
    """
    Newick string with branch lengths. Iterative post-order (safe for deep
    trees). Internal node names included when `with_names`.
    """
    cm = children_map(tree)
    pm = parent_map(tree)

    def bl(nid: str) -> float:
        p = pm.get(nid)
        return float(tree.branch_lengths.get((p, nid), 0.0)) if p is not None else 0.0

    rendered: dict[str, str] = {}
    order: list[str] = []
    stack = [tree.root_id]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(cm.get(n, []))
    for n in reversed(order):
        kids = cm.get(n, [])
        label = _safe(n) if (with_names or not kids) else ""
        if kids:
            inner = ",".join(rendered[c] for c in kids)
            rendered[n] = f"({inner}){label}:{bl(n):.8f}"
        else:
            rendered[n] = f"{label}:{bl(n):.8f}"
    return rendered[tree.root_id] + ";"


def to_ete3(tree: TreeState):
    """Convert to an ete3 Tree (lazy import; for RF/quartet cross-checks)."""
    try:
        from ete3 import Tree
    except ImportError as e:  # pragma: no cover
        raise ImportError("to_ete3 requires ete3 (pip install ete3)") from e
    return Tree(to_newick(tree), format=1)
