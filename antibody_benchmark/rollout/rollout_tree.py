"""Root-conditioned recursive rollout over a fixed observed topology.

Fairness: only root_sequence is conditioned on. Generated parents feed children.
true_sequences must never be passed into this function.
"""

from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING, Mapping, Sequence

from antibody_benchmark.trees import Edge

if TYPE_CHECKING:
    from antibody_benchmark.models.base import EvolutionModel


def rollout_tree(
    model: "EvolutionModel",
    root_sequence: str,
    edges: Sequence[Edge],
    branch_lengths: Mapping[Edge, float],
    seed: int,
    *,
    root_id: str | None = None,
    family_id: str = "",
) -> dict[str, str]:
    """Generate sequences for all nodes from root_sequence + fixed topology/BLs.

    Parameters
    ----------
    model :
        EvolutionModel implementing sample_child(parent_seq, branch_length, rng_seed).
    root_sequence :
        Observed (or naive) root sequence used for conditioning ONLY.
    edges :
        Directed parent→child edges of the observed tree.
    branch_lengths :
        Positive branch length per edge.
    seed :
        Base RNG seed; per-edge seeds are derived deterministically.
    root_id :
        Optional explicit root node id. If omitted, inferred as the unique
        indegree-0 node among ``edges``.
    family_id :
        Optional id mixed into per-edge seeds for reproducibility across families.

    Returns
    -------
    dict[str, str]
        Generated sequence per node id. Does not consult ground-truth non-root
        sequences.
    """
    if not root_sequence:
        raise ValueError("root_sequence must be non-empty")

    children: dict[str, list[str]] = defaultdict(list)
    indeg: dict[str, int] = defaultdict(int)
    nodes: set[str] = set()
    for p, c in edges:
        children[p].append(c)
        indeg[c] += 1
        nodes.add(p)
        nodes.add(c)
        if (p, c) not in branch_lengths:
            raise KeyError(f"missing branch length for edge {(p, c)}")

    if root_id is None:
        roots = [n for n in nodes if indeg[n] == 0]
        if len(roots) != 1:
            raise ValueError(f"cannot infer unique root (nroots={len(roots)})")
        root_id = roots[0]
    elif root_id not in nodes and edges:
        raise ValueError(f"root_id {root_id!r} not present in edges")

    generated: dict[str, str] = {root_id: root_sequence}

    # Preorder / topological from root
    order: list[str] = []
    stack = [root_id]
    while stack:
        u = stack.pop()
        order.append(u)
        for v in reversed(children[u]):
            stack.append(v)

    for parent in order:
        parent_seq = generated[parent]
        for child in children[parent]:
            bl = float(branch_lengths[(parent, child)])
            edge_seed = (
                int(seed)
                + abs(hash((family_id, parent, child))) % 1_000_000
            )
            generated[child] = model.sample_child(
                parent_seq=parent_seq,
                branch_length=bl,
                rng_seed=int(edge_seed % (2**31 - 1)),
            )
    return generated
