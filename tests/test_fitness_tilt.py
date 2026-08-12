"""Unit tests for paper §4.2 fitness tilting (Option A site_local + Option B full_esm)."""

import torch
import torch.nn.functional as F

from src.bridge.fitness_tilt import (
    SCORE_LOG_R0,
    SCORE_LOG_SOFTMAX,
    TILT_FULL_ESM,
    TILT_SITE_LOCAL,
    make_fake_pll_scorer,
    tilt_log_R0_by_fitness,
)
from src.r0_backends import AA_VOCAB


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
    log_R0 = torch.full((1, 2, 20), -2.0)
    log_R0[0, 0, 3] = 2.0
    log_R0[0, 1, 7] = 1.5

    p0 = F.softmax(log_R0, dim=-1)
    p_tilt = F.softmax(
        tilt_log_R0_by_fitness(log_R0, beta=1.0, score=SCORE_LOG_R0), dim=-1
    )

    assert p_tilt[0, 0, 3].item() > p0[0, 0, 3].item()
    assert p_tilt[0, 1, 7].item() > p0[0, 1, 7].item()
    assert p_tilt[0, 0, 0].item() < p0[0, 0, 0].item()


def test_log_softmax_score_mode_also_upweights():
    log_R0 = torch.randn(1, 3, 20)
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


def test_invalid_mode_raises():
    try:
        tilt_log_R0_by_fitness(torch.randn(1, 1, 20), beta=1.0, mode="bogus")
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_full_esm_requires_scorer_and_sequences():
    log_R0 = torch.randn(1, 4, 20)
    try:
        tilt_log_R0_by_fitness(
            log_R0, beta=1.0, mode=TILT_FULL_ESM, sequences=["AAAA"]
        )
        assert False, "expected ValueError without fitness_scorer"
    except ValueError:
        pass


def test_full_esm_tilt_with_fake_scorer():
    """Option B: mutant with high fake-PLL AA should gain mass under β>0."""
    L = 3
    # Base rates nearly uniform.
    log_R0 = torch.zeros(1, L, 20)
    # Fake PLL landscape: AA index 5 is strongly preferred at every site.
    base = torch.full((L, 20), -2.0)
    base[:, 5] = 3.0
    scorer = make_fake_pll_scorer(base)
    seq = "AAA"  # AA_VOCAB index of A is 0
    assert AA_VOCAB[0] == "A"

    p0 = F.softmax(log_R0, dim=-1)
    out = tilt_log_R0_by_fitness(
        log_R0,
        beta=2.0,
        mode=TILT_FULL_ESM,
        sequences=[seq],
        fitness_scorer=scorer,
        batch_size=4,
    )
    p = F.softmax(out, dim=-1)
    # AA 5 (E) should be upweighted vs untilted at each site.
    assert p[0, 0, 5].item() > p0[0, 0, 5].item()
    assert torch.allclose(torch.logsumexp(out, dim=-1), torch.zeros(1, L), atol=1e-5)


def test_full_esm_cache_reused():
    log_R0 = torch.randn(1, 2, 20)
    base = torch.randn(2, 20)
    scorer = make_fake_pll_scorer(base)
    cache: dict = {}
    tilt_log_R0_by_fitness(
        log_R0,
        beta=1.0,
        mode=TILT_FULL_ESM,
        sequences=["AA"],
        fitness_scorer=scorer,
        cache=cache,
        top_k_aas=3,
    )
    n_cached = len(cache)
    assert n_cached > 0
    # Second call should not grow cache for same mutants.
    tilt_log_R0_by_fitness(
        log_R0,
        beta=1.0,
        mode=TILT_FULL_ESM,
        sequences=["AA"],
        fitness_scorer=scorer,
        cache=cache,
        top_k_aas=3,
    )
    assert len(cache) == n_cached


def test_site_local_default_mode():
    log_R0 = torch.randn(1, 2, 20)
    out = tilt_log_R0_by_fitness(log_R0, beta=1.0)  # default mode=site_local
    assert out.shape == log_R0.shape
