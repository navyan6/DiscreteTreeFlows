"""
pLM mutation prior only (row 3): shared BD topology + per-branch ESM mutation.
Not tree-aware — sequences evolve along a birth–death tree using only the
protein language model's per-position distribution.
"""

from __future__ import annotations

import time

from benchmarks.methods.base import Method, GeneratedTree, attach_sequences
from benchmarks.methods.bd_topology import birth_death_topology
from benchmarks.adapters.sequence import evolve_plm


class PLMPrior(Method):
    name = "plm_prior"

    def __init__(self, lm_logits_fn, birth: float, death: float,
                 subs_per_site_scale: float = 1.0):
        """lm_logits_fn(seq) -> [L,20] torch tensor of ESM log-probs."""
        self.lm_logits_fn = lm_logits_fn
        self.birth = birth
        self.death = death
        self.scale = subs_per_site_scale

    def generate(self, root_seq: str, N: int, H: float, seed: int) -> GeneratedTree:
        t0 = time.time()
        topo = birth_death_topology(N, H, self.birth, self.death, seed)
        seqs = evolve_plm(topo, root_seq, self.lm_logits_fn, self.scale, seed)
        tree = attach_sequences(topo, seqs, root_seq)
        return GeneratedTree(tree, {"runtime": time.time() - t0})
