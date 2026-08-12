"""Shared statistical helpers (family-level bootstrap)."""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np


def wasserstein1(a: Sequence[float], b: Sequence[float]) -> float:
    a = np.sort(np.asarray(a, dtype=float))
    b = np.sort(np.asarray(b, dtype=float))
    if a.size == 0 or b.size == 0:
        return float("nan")
    # quantile coupling on a common grid
    qs = np.linspace(0, 1, num=max(len(a), len(b), 50))
    return float(np.mean(np.abs(np.quantile(a, qs) - np.quantile(b, qs))))


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2:
        return float("nan")
    rx = x.argsort().argsort().astype(float)
    ry = y.argsort().argsort().astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = np.sqrt((rx**2).sum() * (ry**2).sum())
    if denom == 0:
        return float("nan")
    return float((rx * ry).sum() / denom)


def js_divergence(p: np.ndarray, q: np.ndarray, eps: float = 1e-12) -> float:
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / (p.sum() + eps)
    q = q / (q.sum() + eps)
    m = 0.5 * (p + q)

    def _kl(a, b):
        mask = a > 0
        return float(np.sum(a[mask] * np.log((a[mask] + eps) / (b[mask] + eps))))

    return 0.5 * _kl(p, m) + 0.5 * _kl(q, m)


def bootstrap_ci(
    values: Sequence[float],
    *,
    n_boot: int = 1000,
    alpha: float = 0.05,
    seed: int = 0,
    agg: Callable[[np.ndarray], float] = np.mean,
) -> dict:
    vals = np.asarray(list(values), dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"point": float("nan"), "lo": float("nan"), "hi": float("nan")}
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        sample = vals[rng.integers(0, vals.size, size=vals.size)]
        boots.append(float(agg(sample)))
    boots_arr = np.sort(boots)
    lo = float(np.quantile(boots_arr, alpha / 2))
    hi = float(np.quantile(boots_arr, 1 - alpha / 2))
    return {"point": float(agg(vals)), "lo": lo, "hi": hi}
