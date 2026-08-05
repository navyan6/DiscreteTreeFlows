"""
Paper §4.2 exponential tilting / fitness weighting (Option A, site-local).

Mutation prior with fitness:
    Q⁰_F(x, x') ∝ Q⁰(x, x') exp(β F(x'))

Cheap site-local proxy (default): treat the ESM site score of amino acid a as
local fitness F for proposing a at that position. With score_a = log R0_a
(raw ESM log-probs / logits) or score_a = log_softmax(log R0)_a:

    log q_F(a) = log_softmax(log_R0)_a + β · score_a
    log R0_tilted = log_softmax(log q_F)     # renormalize

TreeSBM keeps log R_θ = log R0_tilted + c_θ (tilt applied to R0 before RateHeads).

β = 0 disables tilting (identity). Paper-faithful experiments typically use β ≈ 1.
Option B (full mutant PLL) is not implemented — too expensive for train/gen.

Poisson branching λ(x) is intentionally NOT implemented here.
"""

from __future__ import annotations

import torch
import torch.nn.functional as F

SCORE_LOG_R0 = "log_R0"
SCORE_LOG_SOFTMAX = "log_softmax"
_VALID_SCORES = (SCORE_LOG_R0, SCORE_LOG_SOFTMAX)


def tilt_log_R0_by_fitness(
    log_R0: torch.Tensor,
    beta: float = 0.0,
    score: str = SCORE_LOG_R0,
) -> torch.Tensor:
    """
    Apply site-local exponential tilt to ESM reference log-rates.

    Args:
        log_R0: [..., 20] ESM mutation log-rates (typically already log-softmax'd).
        beta: fitness temperature. 0 → return log_R0 unchanged.
        score: fitness proxy per AA —
            "log_R0" (default): score_a = log_R0_a (raw tensor values);
            "log_softmax": score_a = log_softmax(log_R0)_a (normalized site logprobs).

    Returns:
        log_R0_tilted with the same shape, last dim log-normalized when beta ≠ 0.
    """
    if score not in _VALID_SCORES:
        raise ValueError(
            f"fitness score must be one of {_VALID_SCORES}, got {score!r}"
        )
    beta = float(beta)
    if beta == 0.0:
        return log_R0

    # log q(a) = log_softmax(log_R0)_a
    log_q = F.log_softmax(log_R0, dim=-1)
    if score == SCORE_LOG_SOFTMAX:
        score_a = log_q
    else:
        # Default: unnormalized / stored ESM scores (often already logprobs).
        score_a = log_R0

    # log q_F(a) = log q(a) + β · score_a ; then renormalize
    return F.log_softmax(log_q + beta * score_a, dim=-1)
