"""Unified EvolutionModel API for the antibody rollout benchmark."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Mapping


class EvolutionModel(ABC):
    name: str = "base"
    alphabet: str = "nt"  # "nt" or "aa"

    @abstractmethod
    def sample_child(
        self,
        parent_seq: str,
        branch_length: float,
        rng_seed: int,
    ) -> str:
        """Sample a child sequence given parent and branch length."""


class IdentityNullModel(EvolutionModel):
    """TEMP/stub smoke model: returns parent unchanged (loud label in results)."""

    name = "identity_null"
    alphabet = "aa"

    def sample_child(self, parent_seq: str, branch_length: float, rng_seed: int) -> str:
        return parent_seq


class PoissonSiteMutationNullModel(EvolutionModel):
    """TEMP calibration stub: independent AA mutations ~ Poisson(L * t)."""

    name = "poisson_null"
    alphabet = "aa"
    AA = "ACDEFGHIKLMNPQRSTVWY"

    def sample_child(self, parent_seq: str, branch_length: float, rng_seed: int) -> str:
        import numpy as np

        rng = np.random.default_rng(rng_seed)
        seq = list(parent_seq)
        n_mut = int(rng.poisson(max(branch_length, 0.0) * max(len(seq), 1)))
        for _ in range(n_mut):
            if not seq:
                break
            i = int(rng.integers(0, len(seq)))
            choices = [a for a in self.AA if a != seq[i]]
            seq[i] = choices[int(rng.integers(0, len(choices)))]
        return "".join(seq)


def load_models(cfg: Mapping[str, Any]) -> dict[str, EvolutionModel]:
    """Instantiate enabled models from YAML config. Missing deps raise loudly."""
    models: dict[str, EvolutionModel] = {}
    mcfg = cfg.get("models", {})

    if mcfg.get("identity_null", {}).get("enabled", False):
        models["identity_null"] = IdentityNullModel()
        models["poisson_null"] = PoissonSiteMutationNullModel()

    if mcfg.get("thrifty", {}).get("enabled", False):
        from antibody_benchmark.models.thrifty import ThriftyModel

        models["thrifty"] = ThriftyModel(model_name=mcfg["thrifty"].get("name"))

    if mcfg.get("dasm_thrifty", {}).get("enabled", False):
        from antibody_benchmark.models.dasm import DASMThriftyModel

        models["dasm_thrifty"] = DASMThriftyModel(
            dasm_name=mcfg["dasm_thrifty"].get("name"),
            thrifty_name=mcfg["dasm_thrifty"].get("thrifty_name"),
        )

    if mcfg.get("cosine", {}).get("enabled", False):
        from antibody_benchmark.models.cosine import CosineModel

        if mcfg["cosine"].get("guided", False):
            raise RuntimeError(
                "FAIRNESS VIOLATION: CoSiNE Guided Gillespie is forbidden for the "
                "main affinity-maturation benchmark. Set models.cosine.guided=false."
            )
        models["cosine"] = CosineModel(ckpt_path=mcfg["cosine"].get("ckpt"))

    if mcfg.get("treesbm", {}).get("enabled", False):
        from antibody_benchmark.models.treesbm import TreeSBMModel

        models["treesbm"] = TreeSBMModel(
            checkpoint=mcfg["treesbm"].get("checkpoint"),
            force_observed_topology=bool(
                mcfg["treesbm"].get("force_observed_topology", True)
            ),
            max_seq_len=int(mcfg["treesbm"].get("max_seq_len", 566)),
            device=mcfg["treesbm"].get("device"),
            mutation_rate_scale=float(
                mcfg["treesbm"].get("mutation_rate_scale", 1.0)
            ),
        )
    return models
