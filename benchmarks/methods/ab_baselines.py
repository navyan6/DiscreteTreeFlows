"""
Table 5 antibody baselines — wire existing generators; stub true Ab-specific ones.

Paper rows (affinity maturation):
  1. Neutral SHM model          → NOT S5F yet; stub (optional NeutralBD fallback)
  2. pLM mutation prior only    → reuse PLMPrior (ESM-2 + BD topology)
  3. Autoregressive tree edit   → NOT native; stub (ARTreeFormer adapted pools only)

Do not invent S5F / Ab-AR implementations here.
"""

from __future__ import annotations

from benchmarks.methods.base import GeneratedTree, Method
from benchmarks.methods.bd_methods import NeutralBD
from benchmarks.methods.plm_prior import PLMPrior


class NeutralSHMStub(Method):
    """
    Placeholder for Neutral SHM (S5F / hotspot SHM).

    If allow_neutral_bd_fallback=True, delegates to NeutralBD (AA CTMC) so the
    eval harness can run end-to-end — label results TEMP / not paper-final.
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
                "Neutral SHM (S5F) not implemented. Install/wire S5F or pass "
                "allow_neutral_bd_fallback=True for TEMP NeutralBD proxy."
            )
        out = self._fallback.generate(root_seq, N, H, seed)
        out.meta = {**out.meta, "baseline": "neutral_shm_TEMP_neutral_bd_proxy", "paper_ok": False}
        return out


class PLMMutationPriorAb(PLMPrior):
    """Same as Table 2 PLMPrior; renamed for Table 5 row clarity."""

    name = "plm_prior_ab"


class AutoregressiveTreeEditStub(Method):
    """
    Autoregressive tree-edit model stub.

    Native forward Ab AR not in repo. ARTreeFormer adapted pools (viral) are
    not valid Ab baselines — do not silently reuse.
    """

    name = "ar_tree_edit_stub"

    def generate(self, root_seq: str, N: int, H: float, seed: int) -> GeneratedTree:
        raise NotImplementedError(
            "Autoregressive tree-edit baseline for Abs not wired. "
            "Need Ab-trained ARTreeFormer / tree-edit model — see TABLE5_AB_BASELINE_PLAN.md."
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
