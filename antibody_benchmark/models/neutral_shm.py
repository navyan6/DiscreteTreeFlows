"""Simple Neutral SHM baseline (non-Thrifty).

JC69 nucleotide CTMC along the observed topology + branch lengths.
No context-dependent rates (not Thrifty / not S5F neural), no affinity selection.
"""

from __future__ import annotations

import math

from antibody_benchmark.models.base import EvolutionModel


class NeutralSHMModel(EvolutionModel):
    """Independent-site Jukes–Cantor NT evolution (neutral SHM proxy)."""

    name = "neutral_shm"
    alphabet = "nt"
    BASES = "ACGT"

    def sample_child(self, parent_seq: str, branch_length: float, rng_seed: int) -> str:
        import numpy as np

        rng = np.random.default_rng(int(rng_seed) % (2**31 - 1))
        bl = max(float(branch_length), 0.0)
        # JC69: P(site differs) = (3/4) * (1 - exp(-4/3 * t))
        p_change = 0.75 * (1.0 - math.exp(-(4.0 / 3.0) * bl))
        out: list[str] = []
        for b in parent_seq:
            bu = b.upper()
            if bu not in self.BASES:
                out.append(b)
                continue
            if rng.random() < p_change:
                choices = [x for x in self.BASES if x != bu]
                out.append(choices[int(rng.integers(0, 3))])
            else:
                out.append(bu)
        return "".join(out)
