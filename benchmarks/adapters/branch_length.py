"""
Shared branch-length adapter for topology-only methods.

Fit on TRAIN trees only. Assigns branch lengths to a bare topology by sampling
from the empirical pendant/internal branch-length distributions, then rescales
so the mean root-to-tip distance equals the requested H. Same fitted adapter is
used for every topology-only method in the main comparison.
"""

from __future__ import annotations

import json
import random
from pathlib import Path

from src.tree_state import TreeState
from benchmarks.metrics import trees as T


class BranchLengthAdapter:
    def __init__(self):
        self.pendant: list[float] = []
        self.internal: list[float] = []

    def fit(self, train_trees: list[TreeState]) -> "BranchLengthAdapter":
        for t in train_trees:
            leaves = set(T.leaf_labels(t))
            for (p, c), v in t.branch_lengths.items():
                (self.pendant if c in leaves else self.internal).append(float(v))
        return self

    def assign(self, topology: TreeState, H: float, seed: int = 0) -> TreeState:
        rng = random.Random(seed)
        leaves = set(T.leaf_labels(topology))
        pend = self.pendant or [0.01]
        intr = self.internal or [0.01]
        bls = {(p, c): rng.choice(pend if c in leaves else intr)
               for (p, c) in topology.edges}
        ts = TreeState(node_ids=topology.node_ids, root_id=topology.root_id,
                       edges=topology.edges, branch_lengths=bls,
                       node_seqs=topology.node_seqs, active_leaves=list(leaves))
        # rescale to mean root-to-tip == H
        times = T.node_times(ts)
        cur = sum(times[l] for l in leaves) / len(leaves) if leaves else 0.0
        if cur > 0 and H > 0:
            s = H / cur
            bls = {e: v * s for e, v in bls.items()}
            ts = TreeState(node_ids=ts.node_ids, root_id=ts.root_id, edges=ts.edges,
                           branch_lengths=bls, node_seqs=ts.node_seqs,
                           active_leaves=list(leaves))
        return ts

    def save(self, path):
        Path(path).write_text(json.dumps({"pendant": self.pendant, "internal": self.internal}))

    @classmethod
    def load(cls, path):
        d = json.loads(Path(path).read_text())
        a = cls(); a.pendant = d["pendant"]; a.internal = d["internal"]; return a
