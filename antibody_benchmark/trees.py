"""AntibodyTree dataclass and strict validation."""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


Edge = Tuple[str, str]


@dataclass
class AntibodyTree:
    family_id: str
    root_id: str
    node_ids: list[str]
    edges: list[Edge]
    branch_lengths: dict[Edge, float]
    true_sequences: dict[str, str]
    is_leaf: dict[str, bool]
    donor_id: str = ""
    true_aa_sequences: dict[str, str] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    def children_of(self, node_id: str) -> list[str]:
        return [c for p, c in self.edges if p == node_id]

    def preorder(self) -> list[str]:
        order: list[str] = []
        stack = [self.root_id]
        children = defaultdict(list)
        for p, c in self.edges:
            children[p].append(c)
        while stack:
            u = stack.pop()
            order.append(u)
            for v in reversed(children[u]):
                stack.append(v)
        return order

    def leaves(self) -> list[str]:
        return [n for n in self.node_ids if self.is_leaf.get(n, False)]

    @property
    def leaf_ids(self) -> list[str]:
        return self.leaves()

    @property
    def root_sequence(self) -> str:
        return self.true_sequences[self.root_id]


def validate_tree(tree: AntibodyTree) -> tuple[bool, str]:
    """Strict topology/sequence checks for a frozen AntibodyTree.

    Failure reasons (stable strings used in attrition / tests):
      - root_not_in_nodes
      - duplicate_nodes
      - missing_bl / nonpos_bl / invalid_bl
      - edge_unknown_node
      - multi_root / zero_root / root_mismatch
      - multi_parent
      - cycle
      - disconnected
      - missing_seq
      - len_mismatch
    """
    nodes = set(tree.node_ids)
    if tree.root_id not in nodes:
        return False, "root_not_in_nodes"
    if len(tree.node_ids) != len(nodes):
        return False, "duplicate_nodes"

    indeg: Dict[str, int] = defaultdict(int)
    children: Dict[str, List[str]] = defaultdict(list)
    for e in tree.edges:
        if e not in tree.branch_lengths:
            return False, "missing_bl"
        bl = tree.branch_lengths[e]
        if bl is None or isinstance(bl, bool) or not isinstance(bl, (int, float)):
            return False, "invalid_bl"
        if not (float(bl) > 0.0) or float(bl) != float(bl):  # NaN check
            return False, "nonpos_bl"
        p, c = e
        if p not in nodes or c not in nodes:
            return False, "edge_unknown_node"
        indeg[c] += 1
        children[p].append(c)

    roots = [n for n in nodes if indeg[n] == 0]
    if len(roots) == 0:
        return False, "zero_root"
    if len(roots) > 1:
        return False, "multi_root"
    if roots[0] != tree.root_id:
        return False, "root_mismatch"

    # Connectivity from declared root (before cycle/indegree detail)
    seen: set[str] = set()
    q = deque([tree.root_id])
    while q:
        u = q.popleft()
        if u in seen:
            continue
        seen.add(u)
        q.extend(children[u])
    if seen != nodes:
        return False, "disconnected"

    # Directed cycle detection on the connected component
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {n: WHITE for n in nodes}

    def _dfs(u: str) -> bool:
        color[u] = GRAY
        for v in children[u]:
            if color[v] == GRAY:
                return True
            if color[v] == WHITE and _dfs(v):
                return True
        color[u] = BLACK
        return False

    for n in nodes:
        if color[n] == WHITE and _dfs(n):
            return False, "cycle"

    if any(indeg[n] != 1 for n in nodes if n != tree.root_id):
        return False, "multi_parent"

    for p, c in tree.edges:
        sp = tree.true_sequences.get(p, "")
        sc = tree.true_sequences.get(c, "")
        if not sp or not sc:
            return False, "missing_seq"
        if len(sp) != len(sc):
            return False, "len_mismatch"

    return True, "ok"


# Back-compat alias
validate_antibody_tree = validate_tree


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError(f"length mismatch {len(a)} vs {len(b)}")
    return sum(x != y for x, y in zip(a, b))
