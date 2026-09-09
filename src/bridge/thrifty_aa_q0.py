"""
Recipe B: Thrifty (NT SHM) → AA-marginal Q0 for TreeSBM.

OAS clones are amino-acid only. Reverse-translate with human codon usage,
run netam Thrifty on the NT string, collapse codon destination probabilities
to 20 AAs. Output is the same [N, L, 20] log-softmax slot as ESM R0.

Does not overwrite ESM ``group_*_ref_rates.pt`` (cache tag ``_thrifty``).
"""

from __future__ import annotations

from itertools import product
from typing import Optional, Sequence

import torch

from antibody_benchmark.codon import CODON_TABLE
from src.r0_backends import AA_TO_IDX, R0Backend

# Kazusa Homo sapiens codon usage (frequency per 1000). Normalized per AA below.
_HUMAN_CODON_PER_THOUSAND = {
    "TTT": 17.6, "TTC": 20.3, "TTA": 7.7, "TTG": 12.9,
    "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4,
    "TAT": 12.2, "TAC": 15.3, "TAA": 1.0, "TAG": 0.8,
    "TGT": 10.6, "TGC": 12.6, "TGA": 1.6, "TGG": 13.2,
    "CTT": 13.2, "CTC": 19.6, "CTA": 7.2, "CTG": 39.6,
    "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
    "CAT": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2,
    "CGT": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4,
    "ATT": 16.0, "ATC": 20.8, "ATA": 7.5, "ATG": 22.0,
    "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1,
    "AAT": 17.0, "AAC": 19.1, "AAA": 24.4, "AAG": 31.9,
    "AGT": 12.1, "AGC": 19.5, "AGA": 12.2, "AGG": 12.0,
    "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.1,
    "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    "GAT": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6,
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
}

CODON64 = ["".join(p) for p in product("ACGT", repeat=3)]


def _codon_usage_by_aa() -> dict[str, list[tuple[str, float]]]:
    by: dict[str, list[tuple[str, float]]] = {}
    for codon, f in _HUMAN_CODON_PER_THOUSAND.items():
        aa = CODON_TABLE.get(codon, "*")
        if aa == "*":
            continue
        by.setdefault(aa, []).append((codon, float(f)))
    out: dict[str, list[tuple[str, float]]] = {}
    for aa, pairs in by.items():
        tot = sum(p[1] for p in pairs) or 1.0
        out[aa] = [(c, w / tot) for c, w in pairs]
    return out


_USAGE = _codon_usage_by_aa()


def preferred_codon(aa: str) -> str:
    aa = aa.upper()
    if aa not in _USAGE:
        return "NNN"
    return max(_USAGE[aa], key=lambda t: t[1])[0]


def reverse_translate_usage(aa: str) -> str:
    """Deterministic human-usage reverse-translate (argmax codon per AA)."""
    return "".join(preferred_codon(a) if a in _USAGE else "NNN" for a in aa.upper())


def codon_probs_to_aa_logprobs(codon_probs: torch.Tensor) -> torch.Tensor:
    """
    Collapse dest-codon probabilities to 20 AA log-probs.

    Args:
        codon_probs: [n_codons, 64] or [64] (one site).
    Returns:
        [n_codons, 20] log-probs over AA (stops dropped, then renormalized).
    """
    x = codon_probs
    if x.ndim == 1:
        x = x.unsqueeze(0)
    n, k = x.shape
    if k != 64:
        raise ValueError(f"expected 64 codon dests, got {k}")
    aa_mass = x.new_zeros(n, 20)
    for j, codon in enumerate(CODON64):
        aa = CODON_TABLE.get(codon, "*")
        if aa == "*" or aa not in AA_TO_IDX:
            continue
        aa_mass[:, AA_TO_IDX[aa]] = aa_mass[:, AA_TO_IDX[aa]] + x[:, j]
    aa_mass = aa_mass.clamp(min=1e-12)
    aa_mass = aa_mass / aa_mass.sum(dim=-1, keepdim=True)
    return aa_mass.log()


class ThriftyAAR0Backend(R0Backend):
    """
    Frozen Q0 = Thrifty SHM on reverse-translated AA, marginalized to AA.

    ``branch_length`` is the Thrifty CTMC time used to form codon dest probs
    (not TreeSBM sampling BL). Default 0.1 ≈ modest SHM load.
    """

    name = "thrifty_aa"

    def __init__(
        self,
        model_name: str = "ThriftyHumV0.2-59",
        device: Optional[str] = None,
        branch_length: float = 0.1,
    ):
        self.model_name = model_name
        self.branch_length = float(branch_length)
        self.device = device
        self._crepe = None
        self._multihit = None
        self._load()

    def _load(self) -> None:
        try:
            from netam import pretrained
        except ImportError as e:
            raise RuntimeError(
                "Recipe B needs netam (matsengrp/netam). "
                "Install before --r0-backend thrifty_aa."
            ) from e
        self._crepe = pretrained.load(self.model_name, device=self.device)
        try:
            self._multihit = pretrained.load_multihit(
                getattr(self._crepe.model, "multihit_model_name", None)
            )
        except Exception:
            self._multihit = None

    def _codon_probs_nt(self, nt: str) -> torch.Tensor:
        from netam.framework import trimmed_shm_outputs_of_parent_pair
        from netam.molevol import neutral_codon_probs_of_seq
        from netam.sequences import codon_mask_tensor_of

        parent_pair = (nt, "")
        rates, csps = trimmed_shm_outputs_of_parent_pair(self._crepe, parent_pair)
        mask = codon_mask_tensor_of(nt)
        codon_probs = neutral_codon_probs_of_seq(
            nt,
            mask,
            rates[0],
            csps[0],
            self.branch_length,
            multihit_model=self._multihit,
        )
        if not torch.is_tensor(codon_probs):
            codon_probs = torch.as_tensor(codon_probs, dtype=torch.float32)
        else:
            codon_probs = codon_probs.detach().float().cpu()
        if codon_probs.ndim == 1:
            codon_probs = codon_probs.unsqueeze(0)
        if codon_probs.size(-1) != 64:
            raise RuntimeError(
                f"Thrifty codon_probs last dim {tuple(codon_probs.shape)} "
                "expected 64 dest codons"
            )
        return codon_probs

    def log_mutation_rates(
        self,
        sequences: Sequence[str],
        max_seq_len: int,
        device: Optional[torch.device] = None,
    ) -> torch.Tensor:
        N, L = len(sequences), int(max_seq_len)
        log_rates = torch.full((N, L, 20), -2.995732)  # log(1/20)
        for i, seq in enumerate(sequences):
            aa = (seq or "").upper().replace("-", "X")[:L]
            if not aa:
                continue
            nt_chars = set(aa) - set("N.-X")
            if nt_chars and nt_chars <= set("ACGTU"):
                nt = aa.replace("U", "T")
                nt = nt[: len(nt) - len(nt) % 3]
                n_aa = len(nt) // 3
            else:
                clean = "".join(a if a in AA_TO_IDX else "A" for a in aa)
                nt = reverse_translate_usage(clean)
                n_aa = min(len(aa), len(nt) // 3)
            if n_aa <= 0:
                continue
            try:
                cp = self._codon_probs_nt(nt[: n_aa * 3])
            except Exception:
                continue
            aa_lp = codon_probs_to_aa_logprobs(cp[:n_aa])
            clip = min(n_aa, L, aa_lp.size(0))
            log_rates[i, :clip, :] = aa_lp[:clip]
        if device is not None:
            log_rates = log_rates.to(device)
        return log_rates
