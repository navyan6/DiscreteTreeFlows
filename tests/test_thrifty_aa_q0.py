"""Recipe B Q0 helpers (no netam required)."""

from __future__ import annotations

import torch

from src.bridge.thrifty_aa_q0 import (
    CODON64,
    codon_probs_to_aa_logprobs,
    preferred_codon,
    reverse_translate_usage,
)
from src.r0_backends import (
    AA_TO_IDX,
    BACKEND_THRIFTY_AA,
    cache_tag_for_backend,
    list_backends,
    normalize_backend_name,
)


def test_thrifty_aa_registered():
    assert normalize_backend_name("thrifty_aa") == BACKEND_THRIFTY_AA
    assert normalize_backend_name("thrifty_q0") == BACKEND_THRIFTY_AA
    assert cache_tag_for_backend("thrifty_aa") == "_thrifty"
    assert "thrifty_aa" in list_backends()


def test_human_usage_argmax_leu_is_ctg():
    assert preferred_codon("L") == "CTG"
    assert reverse_translate_usage("LM") == "CTGATG"


def test_codon_probs_to_aa_logprobs_stop_dropped_and_renorm():
    p = torch.zeros(64)
    p[CODON64.index("TTT")] = 0.4  # F
    p[CODON64.index("TAA")] = 0.6  # stop
    lp = codon_probs_to_aa_logprobs(p)
    assert lp.shape == (1, 20)
    probs = lp.exp()
    assert abs(probs.sum().item() - 1.0) < 1e-5
    assert abs(probs[0, AA_TO_IDX["F"]].item() - 1.0) < 1e-5
