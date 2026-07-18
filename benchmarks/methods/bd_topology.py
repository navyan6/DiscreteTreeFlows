"""
Shared birth–death topology generator, conditioned on N and H.

Used by the Neutral-CTMC+BD, JTT/WAG/LG+BD, and pLM-prior baselines. Uses
dendropy's birth–death simulator (an established library — no home-rolled BD),
then rescales branch lengths so the mean root-to-tip distance equals H.

Birth/death rates are fitted on TRAIN trees (see benchmarks/fit_params.py) and
passed in; N and H come from the held-out example.
"""

from __future__ import annotations

import random

from src.tree_state import TreeState
from benchmarks.methods.base import dendropy_to_treestate
from benchmarks.metrics import trees as T


def birth_death_topology(N: int, H: float, birth: float, death: float,
                         seed: int = 0, max_tries: int = 50) -> TreeState:
    """Birth–death tree with exactly N extant tips, scaled to mean root-to-tip = H."""
    import dendropy
    from dendropy.simulate import treesim

    rng = random.Random(seed)
    tree = None
    for _ in range(max_tries):
        try:
            tree = treesim.birth_death_tree(
                birth_rate=max(birth, 1e-6),
                death_rate=max(death, 0.0),
                num_extant_tips=N,
                rng=rng,
            )
            break
        except Exception:
            continue
    if tree is None:
        raise RuntimeError(f"birth_death_tree failed for N={N}")

    ts = dendropy_to_treestate(tree)

    # rescale so mean root-to-tip patristic distance == H
    times = T.node_times(ts)
    leaves = T.leaf_labels(ts)
    cur = sum(times[l] for l in leaves) / len(leaves) if leaves else 0.0
    if cur > 0 and H > 0:
        s = H / cur
        ts = TreeState(
            node_ids=ts.node_ids, root_id=ts.root_id, edges=ts.edges,
            branch_lengths={e: v * s for e, v in ts.branch_lengths.items()},
            node_seqs={}, active_leaves=leaves,
        )
    return ts
