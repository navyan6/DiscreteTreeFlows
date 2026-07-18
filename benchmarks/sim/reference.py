"""
Simulated reference track: for a held-out real root, simulate MANY valid
descendant trees under a known evolutionary process (dendropy BD topology +
pyvolve substitution). These give the reference *distribution* that makes
Tree-KL / Split-KL meaningful.

Regimes: neutral, JTT, WAG, LG (birth/death fitted on train).
"""

from __future__ import annotations

from src.tree_state import TreeState
from benchmarks.methods.base import attach_sequences
from benchmarks.methods.bd_topology import birth_death_topology
from benchmarks.adapters.sequence import evolve_pyvolve

REGIMES = ["neutral", "JTT", "WAG", "LG"]


def simulate_reference(root_seq: str, N: int, H: float, regime: str, M: int,
                       birth: float, death: float, seed: int = 0) -> list[TreeState]:
    """M reference descendant trees for one root under `regime`."""
    refs = []
    for m in range(M):
        s = seed * 100003 + m
        topo = birth_death_topology(N, H, birth, death, s)
        seqs = evolve_pyvolve(topo, root_seq, regime, s)
        refs.append(attach_sequences(topo, seqs, root_seq))
    return refs
