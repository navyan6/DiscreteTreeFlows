"""Unit tests for paper §4.2 site-local fitness tilting (Option A)."""

import torch
import torch.nn.functional as F

from src.bridge.fitness_tilt import (
    SCORE_LOG_R0,
    SCORE_LOG_SOFTMAX,
    tilt_log_R0_by_fitness,
)


def test_beta_zero_identity():
    log_R0 = torch.randn(2, 5, 20)
    out = tilt_log_R0_by_fitness(log_R0, beta=0.0)
    assert out is log_R0 or torch.equal(out, log_R0)


def test_beta_zero_score_modes_identity():
    log_R0 = torch.randn(3, 4, 20)
    for score in (SCORE_LOG_R0, SCORE_LOG_SOFTMAX):
        out = tilt_log_R0_by_fitness(log_R0, beta=0.0, score=score)
        assert torch.equal(out, log_R0)


def test_tilted_renormalized():
    log_R0 = torch.randn(2, 7, 20)
    out = tilt_log_R0_by_fitness(log_R0, beta=1.0, score=SCORE_LOG_R0)
    logsum = torch.logsumexp(out, dim=-1)
    assert torch.allclose(logsum, torch.zeros_like(logsum), atol=1e-5)


def test_beta_positive_upweights_high_score_aas():
    """β>0 must increase probability mass on the highest-scoring AA per site."""
    # Hand-crafted: site 0 favors AA 3 strongly under log_R0.
    log_R0 = torch.full((1, 2, 20), -2.0)
    log_R0[0, 0, 3] = 2.0
    log_R0[0, 1, 7] = 1.5

    p0 = F.softmax(log_R0, dim=-1)
    p_tilt = F.softmax(
        tilt_log_R0_by_fitness(log_R0, beta=1.0, score=SCORE_LOG_R0), dim=-1
    )

    assert p_tilt[0, 0, 3].item() > p0[0, 0, 3].item()
    assert p_tilt[0, 1, 7].item() > p0[0, 1, 7].item()
    # Mass on a low-scoring AA should shrink.
    assert p_tilt[0, 0, 0].item() < p0[0, 0, 0].item()


def test_log_softmax_score_mode_also_upweights():
    log_R0 = torch.randn(1, 3, 20)
    # Make AA 5 clearly best at site 1.
    log_R0[0, 1, :] = -3.0
    log_R0[0, 1, 5] = 3.0

    p0 = F.softmax(log_R0, dim=-1)
    p_tilt = F.softmax(
        tilt_log_R0_by_fitness(log_R0, beta=2.0, score=SCORE_LOG_SOFTMAX),
        dim=-1,
    )
    assert p_tilt[0, 1, 5].item() > p0[0, 1, 5].item()


def test_invalid_score_raises():
    try:
        tilt_log_R0_by_fitness(torch.randn(1, 1, 20), beta=1.0, score="pll")
        assert False, "expected ValueError"
    except ValueError:
        pass
