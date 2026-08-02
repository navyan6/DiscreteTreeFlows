"""
Semigroup regularizer L_semi (§4.5).

Paper:
  L_semi = E[ d(G_{s,t}, G_{r,t} ∘ G_{s,r}) ],  0 ≤ s < r < t ≤ t_max
  L_TreeSBM = L_bridge + λ_semi L_semi

Full tree-rollout composition is too expensive every train step. Default mode is
a cheap **rate-space** surrogate: from a fixed bridge state T_s, query the
network at three interval *durations* Δ_st, Δ_sr, Δ_rt (passed as t_scalar),
compose the short-hop rates by time-weighted averaging, and match to the
direct long-hop rates (MSE on mut logits / branch λ / BL / stop).

Optional ``rollout`` distance (sequence-matched RF + Hamming) is exposed for
rare offline checks via ``tree_composition_distance``.
"""

from __future__ import annotations

import random
from typing import Callable

import torch
import torch.nn.functional as F


def sample_time_triple(
    t_max: float = 0.95,
    rng: random.Random | None = None,
    min_gap: float = 1e-3,
) -> tuple[float, float, float]:
    """Sample 0 ≤ s < r < t ≤ t_max with gaps ≥ min_gap when feasible."""
    draw = rng.random if rng is not None else random.random
    # Three ordered uniforms; pad gaps if they collapse.
    xs = sorted(draw() * t_max for _ in range(3))
    s, r, t = xs[0], xs[1], xs[2]
    if t - s < 3 * min_gap:
        # Degenerate draw near zero width — spread across [0, t_max].
        s = 0.0
        t = max(t_max, 3 * min_gap)
        r = 0.5 * t
    else:
        if r - s < min_gap:
            r = s + min_gap
        if t - r < min_gap:
            t = min(t_max, r + min_gap)
            if t <= r:
                r = max(s + min_gap, t - min_gap)
    return float(s), float(r), float(t)


def compose_rates(
    out_sr: dict[str, torch.Tensor],
    out_rt: dict[str, torch.Tensor],
    s: float,
    r: float,
    t: float,
) -> dict[str, torch.Tensor]:
    """Time-weighted linear composition of two short-hop rate predictions."""
    dt = max(t - s, 1e-8)
    w_sr = (r - s) / dt
    w_rt = (t - r) / dt
    return {
        "log_R_theta_mut": w_sr * out_sr["log_R_theta_mut"] + w_rt * out_rt["log_R_theta_mut"],
        "branching_rate": w_sr * out_sr["branching_rate"] + w_rt * out_rt["branching_rate"],
        "branch_length": w_sr * out_sr["branch_length"] + w_rt * out_rt["branch_length"],
        "stop_prob": w_sr * out_sr["stop_prob"] + w_rt * out_rt["stop_prob"],
    }


def semigroup_rate_loss(
    out_st: dict[str, torch.Tensor],
    out_sr: dict[str, torch.Tensor],
    out_rt: dict[str, torch.Tensor],
    s: float,
    r: float,
    t: float,
    w_mut: float = 1.0,
    w_branch: float = 1.0,
    w_bl: float = 1.0,
    w_stop: float = 1.0,
) -> torch.Tensor:
    """
    d(G_st, G_rt ∘ G_sr) in rate space: MSE between direct and composed rates.
    """
    composed = compose_rates(out_sr, out_rt, s, r, t)
    L = w_mut * F.mse_loss(out_st["log_R_theta_mut"], composed["log_R_theta_mut"])
    L = L + w_branch * F.mse_loss(out_st["branching_rate"], composed["branching_rate"])
    L = L + w_bl * F.mse_loss(out_st["branch_length"], composed["branch_length"])
    L = L + w_stop * F.mse_loss(out_st["stop_prob"], composed["stop_prob"])
    return L


def semigroup_loss_from_predictor(
    rates_at_duration: Callable[[float], dict[str, torch.Tensor]],
    s: float,
    r: float,
    t: float,
    **loss_kwargs,
) -> torch.Tensor:
    """
    Query ``rates_at_duration(delta)`` for Δ∈{t−s, r−s, t−r} and return L_semi.

    ``rates_at_duration`` should run tree_enc + rate_heads with t_scalar=delta
    on a fixed bridge state T_s (same features, only time conditioning changes).
    """
    dt_st = max(t - s, 1e-8)
    dt_sr = max(r - s, 1e-8)
    dt_rt = max(t - r, 1e-8)
    out_st = rates_at_duration(dt_st)
    out_sr = rates_at_duration(dt_sr)
    out_rt = rates_at_duration(dt_rt)
    return semigroup_rate_loss(out_st, out_sr, out_rt, s, r, t, **loss_kwargs)


def tree_composition_distance(tree_direct, tree_composed) -> float:
    """
    Optional rollout-mode distance: sequence-matched RF + mean terminal Hamming.

    Used rarely (not every train step). Imports metrics lazily so train path
    stays light when rollout mode is off.
    """
    from benchmarks.metrics.matched import sequence_matched_rf, terminal_edit_distance

    rf = float(sequence_matched_rf(tree_direct, tree_composed))
    te = terminal_edit_distance(tree_direct, tree_composed)["mean"]
    te = float(te) if te == te else 0.0
    return rf + te
