"""
Common interface for every benchmarked method + dendropy→TreeState conversion.

A Method receives only (root_seq, N, H, seed) and returns a GeneratedTree: a
rooted TreeState with branch lengths and internal + terminal sequences. It never
sees the real descendants.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from src.tree_state import TreeState


@dataclass
class GeneratedTree:
    tree: TreeState                 # topology + branch_lengths + node_seqs (internal+terminal)
    meta: dict = field(default_factory=dict)   # runtime, failure info, etc.


class Method(ABC):
    name: str = "method"

    @abstractmethod
    def generate(self, root_seq: str, N: int, H: float, seed: int) -> GeneratedTree:
        ...


def attach_sequences(topology: TreeState, node_seqs: dict[str, str],
                     root_seq: str) -> TreeState:
    """Attach evolved sequences to a topology; force the root to the supplied seq."""
    seqs = {n: node_seqs.get(n, "") for n in topology.node_ids}
    seqs[topology.root_id] = root_seq

    # Some sequence simulators rename unlabeled internal nodes, so we may only
    # get leaf sequences back under stable labels. Backfill missing internals
    # from the nearest ancestor with a sequence so validity checks don't fail on
    # empty strings; leaf-level benchmark metrics only consume terminal seqs.
    children: dict[str, list[str]] = {}
    for parent, child in topology.edges:
        children.setdefault(parent, []).append(child)
    stack = [topology.root_id]
    while stack:
        parent = stack.pop()
        parent_seq = seqs.get(parent, "")
        for child in children.get(parent, []):
            if not seqs.get(child):
                seqs[child] = parent_seq
            stack.append(child)

    return TreeState(node_ids=topology.node_ids, root_id=topology.root_id,
                     edges=topology.edges, branch_lengths=topology.branch_lengths,
                     node_seqs=seqs, active_leaves=list(topology.active_leaves))


def dendropy_to_treestate(t) -> TreeState:
    """Convert a dendropy Tree (topology + branch lengths) to a TreeState (no seqs)."""
    name_map, counter = {}, [0]

    def nm(node):
        if node in name_map:
            return name_map[node]
        if node.taxon is not None and node.taxon.label:
            n = f"L_{node.taxon.label}".replace(" ", "_")
        else:
            n = f"I{counter[0]:06d}"; counter[0] += 1
        name_map[node] = n
        return n

    node_ids, edges, bls = [], [], {}
    for node in t.preorder_node_iter():
        node_ids.append(nm(node))
        for ch in node.child_nodes():
            edges.append((nm(node), nm(ch)))
            bls[(nm(node), nm(ch))] = float(ch.edge.length or 0.0)
    root = nm(t.seed_node)
    leaves = [nm(l) for l in t.leaf_node_iter()]
    return TreeState(node_ids=node_ids, root_id=root, edges=edges,
                     branch_lengths=bls, node_seqs={}, active_leaves=leaves)
