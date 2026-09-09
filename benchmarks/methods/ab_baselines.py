"""
Table 5 antibody baselines — Neutral SHM / pLM / AR for Rod.82 Track C.

Paper rows (affinity maturation):
  1. Neutral SHM model          → JC69 NT CTMC on observed topo+BL (not Thrifty)
  2. pLM mutation prior only    → ESM-2 on observed topo+BL
  3. Autoregressive tree edit   → ARTreeFormer N16 pool pruned + JTT (adapted)

Legacy NeutralBD / stub paths kept for older eval_ab_maturation harness only.
"""

from __future__ import annotations

from benchmarks.methods.base import GeneratedTree, Method
from benchmarks.methods.bd_methods import NeutralBD
from benchmarks.methods.plm_prior import PLMPrior


class NeutralSHMStub(Method):
    """
    Legacy stub for free-gen NeutralBD fallback.

    Prefer antibody_benchmark.models.neutral_shm.NeutralSHMModel for Rod.82.
    """

    name = "neutral_shm_stub"

    def __init__(self, birth: float, death: float, allow_neutral_bd_fallback: bool = False):
        self.birth = birth
        self.death = death
        self.allow_neutral_bd_fallback = allow_neutral_bd_fallback
        self._fallback = NeutralBD(birth, death) if allow_neutral_bd_fallback else None

    def generate(self, root_seq: str, N: int, H: float, seed: int) -> GeneratedTree:
        if self._fallback is None:
            raise NotImplementedError(
                "Use antibody_benchmark NeutralSHMModel (JC69 on observed topo) "
                "via run_rollouts.py --models neutral_shm. "
                "Or pass allow_neutral_bd_fallback=True for TEMP NeutralBD proxy."
            )
        out = self._fallback.generate(root_seq, N, H, seed)
        out.meta = {**out.meta, "baseline": "neutral_shm_TEMP_neutral_bd_proxy", "paper_ok": False}
        return out


class PLMMutationPriorAb(PLMPrior):
    """Same as Table 2 PLMPrior; renamed for Table 5 row clarity (free BD topo)."""

    name = "plm_prior_ab"


class AutoregressiveTreeEditStub(Method):
    """Legacy stub — prefer antibody_benchmark.models.ar_tree_edit.ARTreeEditModel."""

    name = "ar_tree_edit_stub"

    def generate(self, root_seq: str, N: int, H: float, seed: int) -> GeneratedTree:
        raise NotImplementedError(
            "Use antibody_benchmark ARTreeEditModel via run_rollouts.py "
            "--models ar_tree_edit (ARTreeFormer pool + JTT on Rod.82)."
        )


def build_ab_baseline(name: str, birth: float, death: float, esm_logits=None, **kwargs) -> Method:
    key = name.lower().replace("-", "_")
    if key in {"neutral_shm", "neutral_shm_stub"}:
        return NeutralSHMStub(birth, death, allow_neutral_bd_fallback=kwargs.get("allow_neutral_bd_fallback", False))
    if key in {"plm_prior", "plm_prior_ab", "plm"}:
        if esm_logits is None:
            raise ValueError("plm_prior_ab requires esm_logits")
        return PLMMutationPriorAb(esm_logits, birth, death)
    if key in {"ar_tree_edit", "ar_tree_edit_stub", "autoregressive"}:
        return AutoregressiveTreeEditStub()
    raise KeyError(f"Unknown Ab baseline: {name}")
