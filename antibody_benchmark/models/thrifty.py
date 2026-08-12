"""Thrifty (neutral SHM) adapter via official netam pretrained API."""

from __future__ import annotations

from antibody_benchmark.models.base import EvolutionModel


class ThriftyModel(EvolutionModel):
    """Context-dependent SHM without learned affinity selection."""

    name = "thrifty"
    alphabet = "nt"

    def __init__(self, model_name: str | None = "ThriftyHumV0.2-59", device: str | None = None):
        self.model_name = model_name or "ThriftyHumV0.2-59"
        self.device = device
        self._crepe = None
        self._multihit = None
        self._load()

    def _load(self) -> None:
        try:
            from netam import pretrained
        except ImportError as e:
            raise RuntimeError(
                "BLOCKER: netam is not installed. Install matsengrp/netam "
                "(pip install -e antibody_benchmark/data/raw/repos/netam) before "
                "running Thrifty rollouts."
            ) from e
        try:
            self._crepe = pretrained.load(self.model_name, device=self.device)
            # Optional multihit; None is accepted by neutral_codon_probs_of_seq
            try:
                self._multihit = pretrained.load_multihit(
                    getattr(self._crepe.model, "multihit_model_name", None)
                )
            except Exception:
                self._multihit = None
        except Exception as e:
            raise RuntimeError(
                f"BLOCKER: failed to load Thrifty pretrained model {self.model_name!r} "
                f"via netam.pretrained.load. Underlying error: {e}"
            ) from e

    def sample_child(self, parent_seq: str, branch_length: float, rng_seed: int) -> str:
        if self._crepe is None:
            raise RuntimeError("ThriftyModel not loaded")
        try:
            import torch
            from netam.framework import (
                sample_sequence_from_codon_probs,
                trimmed_shm_outputs_of_parent_pair,
            )
            from netam.molevol import neutral_codon_probs_of_seq
            from netam.sequences import codon_mask_tensor_of
        except Exception as e:
            raise RuntimeError(
                "BLOCKER: netam sampling helpers unavailable. "
                f"Underlying error: {e}"
            ) from e

        torch.manual_seed(int(rng_seed))
        # Heavy-only: pair with empty light (official APIs expect length-2 tuples).
        parent_pair = (parent_seq, "")
        rates, csps = trimmed_shm_outputs_of_parent_pair(self._crepe, parent_pair)
        mask = codon_mask_tensor_of(parent_seq)
        codon_probs = neutral_codon_probs_of_seq(
            parent_seq,
            mask,
            rates[0],
            csps[0],
            float(branch_length),
            multihit_model=self._multihit,
        )
        return str(sample_sequence_from_codon_probs(codon_probs))
