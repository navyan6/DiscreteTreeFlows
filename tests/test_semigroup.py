"""Unit tests for the semigroup rate-composition regularizer."""

import random

import torch

from src.bridge.semigroup import (
    compose_rates,
    sample_time_triple,
    semigroup_rate_loss,
    semigroup_loss_from_predictor,
)


def _fake_rates(n=2, L=4, scale=1.0, offset=0.0):
    return {
        "log_R_theta_mut": scale * torch.randn(n, L, 20) + offset,
        "branching_rate": scale * torch.rand(n) + 0.1 + offset,
        "branch_length": scale * torch.rand(n) + 0.01 + offset,
        "stop_prob": torch.sigmoid(scale * torch.randn(n) + offset),
    }


def test_sample_time_triple_ordered():
    rng = random.Random(0)
    for _ in range(50):
        s, r, t = sample_time_triple(t_max=0.95, rng=rng)
        assert 0.0 <= s < r < t <= 0.95 + 1e-9


def test_compose_identity_when_hops_match_direct():
    """If short-hop rates equal the direct rates, composition reproduces them."""
    out = _fake_rates(scale=0.0, offset=1.0)  # constant rates
    s, r, t = 0.1, 0.4, 0.9
    composed = compose_rates(out, out, s, r, t)
    for k in out:
        assert torch.allclose(composed[k], out[k], atol=1e-6)


def test_infinitesimal_gap_near_zero_loss():
    """When Δt→0 and all three queries return the same rates, L_semi≈0."""
    out = _fake_rates(scale=0.5)
    s, r, t = 0.5, 0.5 + 1e-4, 0.5 + 2e-4
    loss = semigroup_rate_loss(out, out, out, s, r, t)
    assert float(loss) < 1e-6


def test_finite_gap_when_rates_disagree():
    """Disagreeing short-hop vs direct rates yields strictly positive loss."""
    out_st = _fake_rates(scale=1.0, offset=0.0)
    out_sr = _fake_rates(scale=1.0, offset=3.0)
    out_rt = _fake_rates(scale=1.0, offset=-2.0)
    s, r, t = 0.1, 0.4, 0.9
    loss = semigroup_rate_loss(out_st, out_sr, out_rt, s, r, t)
    assert float(loss) > 0.1


def test_predictor_wrapper_matches_manual():
    cache = {
        0.8: _fake_rates(offset=0.0),
        0.3: _fake_rates(offset=1.0),
        0.5: _fake_rates(offset=2.0),
    }

    def rates_at(delta):
        # Exact key lookup for the fixed triple below.
        key = round(delta, 1)
        return cache[key]

    s, r, t = 0.1, 0.4, 0.9  # durations 0.8, 0.3, 0.5
    wrapped = semigroup_loss_from_predictor(rates_at, s, r, t)
    manual = semigroup_rate_loss(cache[0.8], cache[0.3], cache[0.5], s, r, t)
    assert torch.allclose(wrapped, manual)


def test_loss_differentiable():
    out_st = {k: v.clone().requires_grad_(True) for k, v in _fake_rates().items()}
    out_sr = _fake_rates(offset=0.5)
    out_rt = _fake_rates(offset=-0.5)
    loss = semigroup_rate_loss(out_st, out_sr, out_rt, 0.1, 0.4, 0.9)
    loss.backward()
    assert out_st["log_R_theta_mut"].grad is not None
