"""
Native birth–death + substitution baselines (rows 1–2 of the table):
  Neutral CTMC + BD  — pyvolve "neutral" (uniform AA CTMC)
  JTT/WAG/LG + BD    — pyvolve empirical AA models

Same BD topology procedure (bd_topology), differing only in the substitution
model. Birth/death rates fitted on train (passed in).
"""

from __future__ import annotations

import time

from benchmarks.methods.base import Method, GeneratedTree, attach_sequences
from benchmarks.methods.bd_topology import birth_death_topology
from benchmarks.adapters.sequence import evolve_pyvolve


class BDMethod(Method):
    def __init__(self, name: str, model: str, birth: float, death: float):
        self.name = name
        self.model = model
        self.birth = birth
        self.death = death

    def generate(self, root_seq: str, N: int, H: float, seed: int) -> GeneratedTree:
        t0 = time.time()
        topo = birth_death_topology(N, H, self.birth, self.death, seed)
        seqs = evolve_pyvolve(topo, root_seq, self.model, seed)
        tree = attach_sequences(topo, seqs, root_seq)
        return GeneratedTree(tree, {"runtime": time.time() - t0, "model": self.model})


class NeutralBD(BDMethod):
    def __init__(self, birth: float, death: float):
        super().__init__("neutral_bd", "neutral", birth, death)


class EmpiricalBD(BDMethod):
    """JTT/WAG/LG + BD. `model` picked per run or best-on-val for the final table."""
    def __init__(self, birth: float, death: float, model: str = "JTT"):
        super().__init__("empirical_bd", model, birth, death)
