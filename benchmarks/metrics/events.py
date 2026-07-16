"""
Event-head calibration metrics (Track A3 + "calibration of branching/mutation
probabilities").

These score the rate heads directly — predicted next-event distributions and
timings against the true next event read off a real tree — rather than waiting
to judge a whole generated tree. All pure numpy; the driver (track_a) supplies
predictions and true labels.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "event_nll", "brier_score", "expected_calibration_error",
    "time_to_event_error", "mutation_site_recall_at_k",
]


def event_nll(pred_probs: np.ndarray, true_idx: np.ndarray, eps: float = 1e-12) -> float:
    """
    Mean categorical negative log-likelihood of the true event class.
    pred_probs: [N, C] (rows sum to 1); true_idx: [N] in [0, C).
    """
    pred_probs = np.asarray(pred_probs, float)
    true_idx = np.asarray(true_idx, int)
    p = pred_probs[np.arange(len(true_idx)), true_idx]
    return float(-np.mean(np.log(np.clip(p, eps, 1.0))))


def brier_score(pred_probs: np.ndarray, true_idx: np.ndarray) -> float:
    """Multiclass Brier score: mean over samples of sum_c (p_c - 1[c=true])^2."""
    pred_probs = np.asarray(pred_probs, float)
    true_idx = np.asarray(true_idx, int)
    onehot = np.zeros_like(pred_probs)
    onehot[np.arange(len(true_idx)), true_idx] = 1.0
    return float(np.mean(np.sum((pred_probs - onehot) ** 2, axis=1)))


def expected_calibration_error(confidences: np.ndarray, correct: np.ndarray,
                               n_bins: int = 10) -> float:
    """
    ECE over equal-width confidence bins:
      sum_b (n_b / N) * |accuracy_b - mean_confidence_b|.
    confidences: [N] predicted top-class prob; correct: [N] bool (top-class right).
    """
    confidences = np.asarray(confidences, float)
    correct = np.asarray(correct, float)
    N = len(confidences)
    if N == 0:
        return float("nan")
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        # include right edge in the last bin
        in_bin = (confidences > lo) & (confidences <= hi) if hi < 1.0 \
            else (confidences > lo) & (confidences <= hi + 1e-9)
        nb = int(in_bin.sum())
        if nb == 0:
            continue
        acc = correct[in_bin].mean()
        conf = confidences[in_bin].mean()
        ece += (nb / N) * abs(acc - conf)
    return float(ece)


def time_to_event_error(pred_times: np.ndarray, true_times: np.ndarray) -> dict[str, float]:
    """MAE and RMSE between predicted and true waiting times to the next event."""
    pred_times = np.asarray(pred_times, float)
    true_times = np.asarray(true_times, float)
    err = pred_times - true_times
    return {"mae": float(np.mean(np.abs(err))),
            "rmse": float(np.sqrt(np.mean(err ** 2)))}


def mutation_site_recall_at_k(site_scores: np.ndarray, true_sites, ks=(1, 5, 10)
                              ) -> dict[int, float]:
    """
    For one intermediate state: of the sites that actually mutate on the next
    event(s), what fraction fall in the model's top-k highest mutation-propensity
    sites. site_scores: [L] per-site mutation propensity; true_sites: iterable of
    site indices that truly mutate. Returns {k: recall@k}.
    """
    site_scores = np.asarray(site_scores, float)
    true_sites = set(int(s) for s in true_sites)
    out: dict[int, float] = {}
    if not true_sites:
        return {k: float("nan") for k in ks}
    order = np.argsort(-site_scores)  # descending propensity
    for k in ks:
        topk = set(int(i) for i in order[:k])
        out[k] = len(topk & true_sites) / len(true_sites)
    return out
