"""GY94-style codon generator collapsed to a 20×20 AA jump chain (viral R0).

Does not overwrite ESM ``group_*_ref_rates.pt``. Cache tag ``_codon_gy94``.
True CDS is unused here: each AA is a uniform mixture over its sense codons
(standard code). Destinations follow Goldman–Yang (κ transitions, ω nonsyn).
"""

from __future__ import annotations

import numpy as np
import torch

from src.r0_backends import AA_TO_IDX, AA_VOCAB, R0Backend

# Standard genetic code, sense only (no stops).
_CODON_TO_AA: dict[str, str] = {}
_AA_TO_CODONS: dict[str, list[str]] = {aa: [] for aa in AA_VOCAB}

_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y",
    "TGT": "C", "TGC": "C", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
for _c, _a in _TABLE.items():
    if _a in AA_TO_IDX:
        _CODON_TO_AA[_c] = _a
        _AA_TO_CODONS[_a].append(_c)

_NT = "ACGT"
_TS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def _nt_diffs(a: str, b: str) -> list[tuple[str, str]]:
    return [(x, y) for x, y in zip(a, b) if x != y]


def gy94_aa_rate_matrix(kappa: float = 2.0, omega: float = 0.25) -> np.ndarray:
    """Instantaneous 20×20 AA rates by mixing GY94 codon Q (uniform π)."""
    Qaa = np.zeros((20, 20), dtype=np.float64)
    for aa_i, codons_i in _AA_TO_CODONS.items():
        i = AA_TO_IDX[aa_i]
        n_i = max(len(codons_i), 1)
        for c in codons_i:
            for d, aa_j in _CODON_TO_AA.items():
                diffs = _nt_diffs(c, d)
                if len(diffs) != 1:
                    continue
                x, y = diffs[0]
                rate = kappa if (x, y) in _TS else 1.0
                if aa_j != aa_i:
                    rate *= omega
                Qaa[i, AA_TO_IDX[aa_j]] += rate / n_i
    for i in range(20):
        Qaa[i, i] = 0.0
        Qaa[i, i] = -Qaa[i].sum()
    return Qaa


class CodonGY94R0Backend(R0Backend):
    """Sitewise dest log-probs from GY94 collapsed to AA (context-free)."""

    name = "codon_gy94"

    def __init__(self, kappa: float = 2.0, omega: float = 0.25, stay_mass: float = 0.5):
        self.kappa = float(kappa)
        self.omega = float(omega)
        self.stay_mass = float(stay_mass)
        self._Q = gy94_aa_rate_matrix(self.kappa, self.omega)

    def log_mutation_rates(self, sequences, max_seq_len, device=None):
        N, L = len(sequences), max_seq_len
        log_rates = torch.zeros(N, L, 20, dtype=torch.float32)
        stay = min(max(self.stay_mass, 1e-6), 1.0 - 1e-6)
        Q = self._Q
        for i, seq in enumerate(sequences):
            for pos, aa in enumerate(seq[:L]):
                a = AA_TO_IDX.get(aa)
                if a is None:
                    log_rates[i, pos, :] = -np.log(20.0)
                    continue
                row = Q[a].copy()
                off = np.maximum(row, 0.0)
                off[a] = 0.0
                s = off.sum()
                probs = np.full(20, (1.0 - stay) / 19.0, dtype=np.float64)
                if s > 1e-12:
                    probs = off / s * (1.0 - stay)
                probs[a] = stay
                probs = np.clip(probs, 1e-12, None)
                probs /= probs.sum()
                log_rates[i, pos, :] = torch.from_numpy(np.log(probs)).float()
        return log_rates
