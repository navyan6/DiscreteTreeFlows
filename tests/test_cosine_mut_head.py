"""CoSiNE-style mutation head: shapes, no-residual default mix, residual default intact."""

from __future__ import annotations

import torch

from src.bridge.cosine_mut_head import CosineSiteMutHead, expected_init_mix
from src.networks import RateHeads


def test_cosine_site_head_shape_and_init_mix():
    head = CosineSiteMutHead(d_model=8, max_seq_len=16)
    n, L = 2, 10
    h = torch.randn(n, 8)
    aa = torch.randint(0, 20, (n, L))
    log_r0 = torch.randn(n, L, 20)
    out = head(h, aa, log_r0)
    assert out.shape == (n, L, 20)
    mix = float(torch.sigmoid(head.r0_mix_logit))
    assert abs(mix - expected_init_mix(-4.0)) < 1e-5
    assert mix < 0.05


def test_rate_heads_cosine_does_not_equal_r0():
    heads = RateHeads(d_model=16, max_seq_len=12, mut_head_type="cosine")
    H = torch.randn(3, 16)
    log_r0 = torch.zeros(2, 8, 20)
    aa = torch.randint(0, 20, (2, 8))
    out = heads(H, [0, 2], log_r0, aa_indices=aa)
    assert out["log_R_theta_mut"].shape == (2, 8, 20)
    assert not torch.allclose(out["log_R_theta_mut"], log_r0)
    assert heads.needs_aa_indices
    try:
        heads(H, [0, 2], log_r0, aa_indices=None)
        raise AssertionError("expected ValueError")
    except ValueError:
        pass


def test_rate_heads_residual_default_still_adds_c_theta():
    heads = RateHeads(d_model=16, max_seq_len=12)
    assert heads.mut_head_type == "residual"
    assert heads.needs_aa_indices is False
    H = torch.randn(3, 16)
    log_r0 = torch.randn(2, 8, 20)
    out = heads(H, [0, 1], log_r0)
    assert out["log_R_theta_mut"].shape == (2, 8, 20)
