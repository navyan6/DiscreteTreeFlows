"""
Thin adapter: score HA AA sequences with a frozen VaxSeer dominance LM.

TreeSBM remains the generator; VaxSeer is an external scorer only (no retrain).

Setup (Betty / LABHOME):
  git clone https://github.com/wxsh1213/vaxseer.git $LABHOME/third_party/vaxseer
  cd $LABHOME/third_party/vaxseer
  conda env create -f environment.yaml   # or isolate deps
  python download_models_from_dropbox.py --task lm --year <Y> --subtype a_h3n2 --output_dir runs
  # repeat for a_h1n1 year-matched subsets (full dump ~60GB — prefer year subsets)

Env:
  VAXSEER_ROOT  path to cloned repo (default: $LABHOME/third_party/vaxseer)
  VAXSEER_RUNS  path to runs/ with lm weights
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Optional, Sequence


def default_vaxseer_root() -> Path:
    lab = os.environ.get("LABHOME", "")
    if lab:
        return Path(lab) / "third_party" / "vaxseer"
    return Path("third_party") / "vaxseer"


class VaxSeerDominanceScorer:
    """Score AA HA strings → dominance p_t(x) for a season/subtype."""

    def __init__(
        self,
        subtype: str = "a_h3n2",
        year: int = 2024,
        vaxseer_root: Optional[str | Path] = None,
        runs_dir: Optional[str | Path] = None,
        device: str = "cpu",
    ):
        self.subtype = subtype
        self.year = int(year)
        self.device = device
        self.root = Path(vaxseer_root or os.environ.get("VAXSEER_ROOT") or default_vaxseer_root())
        self.runs = Path(
            runs_dir
            or os.environ.get("VAXSEER_RUNS")
            or (self.root / "runs")
        )
        self._model = None
        self._ready = False
        self._error: Optional[str] = None

    def _ensure_path(self) -> None:
        if str(self.root) not in sys.path:
            sys.path.insert(0, str(self.root))

    def load(self) -> None:
        """Best-effort load of official lm weights. Raises with setup hints on failure."""
        self._ensure_path()
        if not self.root.is_dir():
            raise FileNotFoundError(
                f"VaxSeer repo not found at {self.root}. "
                "Clone: git clone https://github.com/wxsh1213/vaxseer.git"
            )
        # Official package layout varies by release; try common entry points.
        try:
            # Prefer a lightweight import of their LM inference if present.
            from vaxseer.dominance import load_dominance_model  # type: ignore

            self._model = load_dominance_model(
                subtype=self.subtype, year=self.year, runs_dir=str(self.runs), device=self.device
            )
            self._ready = True
            return
        except Exception as e1:
            self._error = f"vaxseer.dominance import failed: {e1}"
        try:
            # Fallback: user-facing script API (if they expose score_sequences).
            import importlib

            mod = importlib.import_module("bin.score_dominance")  # type: ignore
            self._model = mod
            self._ready = True
            return
        except Exception as e2:
            self._error = f"{self._error}; bin.score_dominance failed: {e2}"
        raise RuntimeError(
            "Could not load VaxSeer dominance LM. "
            f"root={self.root} runs={self.runs} err={self._error}. "
            "Download: python download_models_from_dropbox.py --task lm "
            f"--year {self.year} --subtype {self.subtype} --output_dir runs"
        )

    def score(self, sequences: Sequence[str], t: Optional[float] = None) -> list[float]:
        """
        Return p_t(x) for each AA HA string.

        If the official API is unavailable, raises RuntimeError (do not invent scores).
        """
        if not self._ready:
            self.load()
        assert self._model is not None
        if hasattr(self._model, "score"):
            return list(self._model.score(list(sequences), t=t))
        if hasattr(self._model, "score_sequences"):
            return list(self._model.score_sequences(list(sequences), year=self.year))
        raise RuntimeError("Loaded VaxSeer object has no score / score_sequences API")


def compare_generated_vs_observed(
    generated: dict[str, Sequence[str]],
    observed: Sequence[str],
    observed_freqs: Optional[Sequence[float]] = None,
    scorer: Optional[VaxSeerDominanceScorer] = None,
    eps_list: Sequence[int] = (2, 5),
) -> dict:
    """
    Dominance recovery + distribution summary.

    generated: method_name -> list of AA leaves
    observed: circulating HA in season Y
    """
    import numpy as np

    if scorer is None:
        raise ValueError("scorer required")
    p_obs = scorer.score(observed)
    p_obs_arr = np.asarray(p_obs, dtype=float)
    thr = float(np.quantile(p_obs_arr, 0.9)) if len(p_obs_arr) else 0.0
    top_idx = [i for i, p in enumerate(p_obs_arr) if p >= thr]
    top_seqs = [observed[i] for i in top_idx]

    def hamming(a: str, b: str) -> int:
        L = min(len(a), len(b))
        return sum(1 for i in range(L) if a[i] != b[i]) + abs(len(a) - len(b))

    out: dict = {
        "n_observed": len(observed),
        "n_top_decile": len(top_seqs),
        "p_obs_mean": float(p_obs_arr.mean()) if len(p_obs_arr) else None,
        "p_obs_median": float(np.median(p_obs_arr)) if len(p_obs_arr) else None,
        "p_obs_max": float(p_obs_arr.max()) if len(p_obs_arr) else None,
        "methods": {},
    }
    for method, gens in generated.items():
        p_g = scorer.score(list(gens))
        p_g_arr = np.asarray(p_g, dtype=float)
        rec = {}
        for eps in eps_list:
            hit = 0
            for o in top_seqs:
                if any(hamming(g, o) <= eps for g in gens):
                    hit += 1
            rec[f"top_decile_recovery_eps{eps}"] = (
                hit / len(top_seqs) if top_seqs else None
            )
        out["methods"][method] = {
            "n_generated": len(gens),
            "p_mean": float(p_g_arr.mean()) if len(p_g_arr) else None,
            "p_median": float(np.median(p_g_arr)) if len(p_g_arr) else None,
            "p_max": float(p_g_arr.max()) if len(p_g_arr) else None,
            **rec,
        }
    if observed_freqs is not None:
        out["observed_freq_provided"] = True
    return out


def load_precomputed_lm_csv(csv_path: str | Path) -> dict[str, float]:
    """Map GISAID EPI id -> dominance score from official results dump.

    Example path under LABHOME/third_party/vaxseer_results:
      runs/pipeline/2018-02/a_h3n2/.../dominance_prediction/lm/.../test_results.csv
    Does **not** score novel AA strings (need live lm weights for that).
    """
    import csv
    out: dict[str, float] = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sid = row.get("src_id") or row.get("id") or ""
            if not sid:
                continue
            try:
                out[sid] = float(row.get("prediction") or row.get("prob") or 0.0)
            except ValueError:
                continue
    return out
