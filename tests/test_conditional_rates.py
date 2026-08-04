"""Unit tests for conditional bridge rate targets (Doob-tilt surrogate)."""

import torch
import torch.nn.functional as F

from src.bridge.conditional_rates import (
    conditional_bridge_kl,
    conditional_bridge_log_target,
)


def test_target_normalized_mid_t():
    torch.manual_seed(0)
    log_R0 = torch.randn(3, 7, 20)
    x1 = torch.randint(0, 20, (3, 7))
    lt = conditional_bridge_log_target(log_R0, x1, t=0.5, c=1.0)
    p = lt.exp()
    assert torch.allclose(p.sum(-1), torch.ones(3, 7), atol=1e-5)
    assert (p >= 0).all()
    assert torch.isfinite(lt).all()


def test_near_one_and_exact_one_finite():
    """P0: t=0.9999 and t=1.0 must not NaN; t=1 is pure CE on x1."""
    torch.manual_seed(1)
    n, L = 2, 5
    log_R0 = torch.randn(n, L, 20)
    x1 = torch.randint(0, 20, (n, L))
    log_theta = torch.randn(n, L, 20)

    for t in (0.9999, 1.0):
        lt = conditional_bridge_log_target(log_R0, x1, t=t, c=1.0)
        assert torch.isfinite(lt).all(), f"non-finite log_target at t={t}"
        p = lt.exp()
        assert torch.isfinite(p).all(), f"non-finite target mass at t={t}"
        assert torch.allclose(p.sum(-1), torch.ones(n, L), atol=1e-5)
        mass_on_x1 = p.gather(-1, x1.unsqueeze(-1)).squeeze(-1)
        assert mass_on_x1.min() > 0.99, f"should concentrate on x1 at t={t}"

        kl = conditional_bridge_kl(log_theta, log_R0, x1, t=t, c=1.0)
        assert torch.isfinite(kl).all(), f"non-finite KL at t={t}"

    # Exact t=1: KL equals per-site CE
    kl1 = conditional_bridge_kl(log_theta, log_R0, x1, t=1.0, c=1.0)
    ce1 = F.cross_entropy(
        log_theta.reshape(-1, 20), x1.reshape(-1), reduction="none"
    ).reshape(n, L)
    assert torch.allclose(kl1, ce1, atol=1e-4)


def test_t_greater_than_one_treated_as_terminal():
    torch.manual_seed(2)
    log_R0 = torch.randn(1, 3, 20)
    x1 = torch.randint(0, 20, (1, 3))
    lt = conditional_bridge_log_target(log_R0, x1, t=1.5, c=1.0)
    assert torch.isfinite(lt).all()
    mass = lt.exp().gather(-1, x1.unsqueeze(-1)).squeeze(-1)
    assert torch.allclose(mass, torch.ones_like(mass), atol=1e-5)


def test_kl_approaches_ce_as_t_to_one():
    torch.manual_seed(3)
    n, L = 2, 4
    log_R0 = torch.randn(n, L, 20)
    x1 = torch.randint(0, 20, (n, L))
    log_theta = torch.randn(n, L, 20)
    gaps = []
    for t in (0.9, 0.99, 0.999, 0.9999):
        kl = conditional_bridge_kl(log_theta, log_R0, x1, t=t, c=1.0)
        ce = F.cross_entropy(
            log_theta.reshape(-1, 20), x1.reshape(-1), reduction="none"
        ).reshape(n, L)
        gaps.append((kl - ce).abs().max().item())
    assert gaps == sorted(gaps, reverse=True)
    assert gaps[-1] < 2e-3
