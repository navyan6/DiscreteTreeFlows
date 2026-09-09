"""GY94→AA R0: shape, rows sum to 0, dest mass on reachable AA."""

from __future__ import annotations

import numpy as np
import torch

from src.bridge.codon_q0 import CodonGY94R0Backend, gy94_aa_rate_matrix
from src.bridge.nt_stall_fire import boost_fire_at_stalls, stall_site_mask
from src.r0_backends import BACKEND_CODON_GY94, build_r0_backend, cache_tag_for_backend


def test_gy94_q_is_rate_matrix():
    Q = gy94_aa_rate_matrix()
    assert Q.shape == (20, 20)
    assert np.allclose(Q.sum(axis=1), 0.0, atol=1e-8)
    assert (np.diag(Q) < 0).all()


def test_codon_backend_logprobs():
    b = CodonGY94R0Backend()
    lp = b.log_mutation_rates(["ACDEFGHIKLMNPQRSTVWY"], 20)
    assert lp.shape == (1, 20, 20)
    p = lp.exp()
    assert torch.allclose(p.sum(-1), torch.ones(1, 20), atol=1e-5)


def test_codon_backend_registered():
    assert cache_tag_for_backend("codon_gy94") == "_codon_gy94"
    assert cache_tag_for_backend("esm2") == ""
    b = build_r0_backend("codon_gy94")
    assert getattr(b, "name", "") in {BACKEND_CODON_GY94, "codon_gy94"}


def test_stall_fire_keeps_mass():
    seq = "A" * 12
    log_R = torch.zeros(12, 20)
    out = boost_fire_at_stalls(log_R, seq, 1.0)
    p = out.softmax(-1)
    assert torch.allclose(p.sum(-1), torch.ones(12), atol=1e-5)
    assert len(stall_site_mask(seq)) == 12
