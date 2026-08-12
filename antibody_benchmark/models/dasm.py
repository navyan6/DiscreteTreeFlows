"""DASM + Thrifty adapter via official netam APIs.

Named DASM+Thrifty because generation combines Thrifty neutral SHM rates with
DASM amino-acid selection factors (Matsen et al.).
"""

from __future__ import annotations

from antibody_benchmark.models.base import EvolutionModel


class DASMThriftyModel(EvolutionModel):
    name = "dasm_thrifty"
    alphabet = "nt"

    def __init__(
        self,
        dasm_name: str | None = "DASMHumV1.0-4M",
        thrifty_name: str | None = "ThriftyHumV0.2-59",
        device: str | None = None,
    ):
        self.dasm_name = dasm_name or "DASMHumV1.0-4M"
        self.thrifty_name = thrifty_name or "ThriftyHumV0.2-59"
        self.device = device
        self._dasm = None
        self._thrifty = None
        self._multihit = None
        self._load()

    def _load(self) -> None:
        try:
            from netam import pretrained
        except ImportError as e:
            raise RuntimeError(
                "BLOCKER: netam is not installed. Install matsengrp/netam before "
                "running DASM+Thrifty rollouts."
            ) from e
        try:
            self._dasm = pretrained.load(self.dasm_name, device=self.device)
            self._thrifty = pretrained.load(self.thrifty_name, device=self.device)
            try:
                self._multihit = pretrained.load_multihit(
                    getattr(self._dasm.model, "multihit_model_name", None)
                )
            except Exception:
                self._multihit = None
        except Exception as e:
            raise RuntimeError(
                f"BLOCKER: failed to load DASM={self.dasm_name!r} / "
                f"Thrifty={self.thrifty_name!r} via netam.pretrained. Error: {e}"
            ) from e

    def sample_child(self, parent_seq: str, branch_length: float, rng_seed: int) -> str:
        if self._dasm is None or self._thrifty is None:
            raise RuntimeError("DASMThriftyModel not loaded")
        try:
            import torch
            from netam.framework import (
                codon_probs_of_parent_seq,
                sample_sequence_from_codon_probs,
            )
        except Exception as e:
            raise RuntimeError(
                f"BLOCKER: netam.framework helpers unavailable: {e}"
            ) from e

        torch.manual_seed(int(rng_seed))
        parent_pair = (parent_seq, "")
        heavy_codon_probs, _light = codon_probs_of_parent_seq(
            self._dasm,
            parent_pair,
            float(branch_length),
            neutral_crepe=self._thrifty,
            multihit_model=self._multihit,
        )
        return str(sample_sequence_from_codon_probs(heavy_codon_probs))
