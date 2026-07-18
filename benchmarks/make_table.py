#!/usr/bin/env python3
"""
Aggregate the long-format benchmark results into the paper table(s).

Input: a long-format CSV where each row is one (method, track, N, H_bucket,
root_id) with the six metric columns already reduced over the K samples for that
root. Output: per-(track, condition) CSV + LaTeX tables with mean ± SE (or
bootstrap 95% CI) over roots, and a main table averaged across conditions.

    python benchmarks/make_table.py --results benchmarks/results/results.csv \
        --out benchmarks/results/tables --ci se
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = ["tree_kl", "split_kl", "rf", "quartet", "branch_w_all", "terminal_edit"]
METRIC_LABELS = {
    "tree_kl": "Tree-KL", "split_kl": "Split-KL", "rf": "RF",
    "quartet": "Quartet", "branch_w_all": "Branch-W1", "terminal_edit": "Term-edit",
}
METHOD_ORDER = [
    "neutral_bd", "empirical_bd", "plm_prior", "artreeformer_adapted",
    "phylovae_adapted", "treesbm",
]


def _agg(vals: np.ndarray, ci: str) -> tuple[float, float]:
    vals = vals[~np.isnan(vals)]
    if len(vals) == 0:
        return float("nan"), float("nan")
    m = float(vals.mean())
    if ci == "boot" and len(vals) > 1:
        boot = [np.random.choice(vals, len(vals), replace=True).mean() for _ in range(1000)]
        lo, hi = np.percentile(boot, [2.5, 97.5])
        return m, float((hi - lo) / 2)     # half-width
    se = float(vals.std(ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
    return m, se


def summarize(df: pd.DataFrame, ci: str) -> pd.DataFrame:
    rows = []
    for method in [m for m in METHOD_ORDER if m in set(df["method"])]:
        sub = df[df["method"] == method]
        row = {"method": method, "n_roots": sub["root_id"].nunique()}
        for met in METRICS:
            if met in sub:
                m, u = _agg(sub[met].to_numpy(dtype=float), ci)
                row[met] = m
                row[met + "_unc"] = u
        rows.append(row)
    return pd.DataFrame(rows)


def to_latex(summary: pd.DataFrame, caption: str) -> str:
    header = " & ".join(["Method"] + [METRIC_LABELS[m] + " $\\downarrow$" for m in METRICS])
    lines = [r"\begin{table}[t]\centering", f"\\caption{{{caption}}}",
             r"\begin{tabular}{l" + "c" * len(METRICS) + "}", r"\toprule", header + r" \\", r"\midrule"]
    for _, r in summary.iterrows():
        cells = [r["method"].replace("_", "\\_")]
        for m in METRICS:
            if m in r and r[m] == r[m]:
                cells.append(f"{r[m]:.3f} $\\pm$ {r[m + '_unc']:.3f}")
            else:
                cells.append("--")
        lines.append(" & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", default="benchmarks/results/tables")
    ap.add_argument("--ci", choices=["se", "boot"], default="se")
    args = ap.parse_args()

    df = pd.read_csv(args.results)
    if "valid" in df:
        df = df[df["valid"] == True]  # noqa: E712
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    # per-(track, N) tables + a main table averaged across conditions per track
    for track in sorted(df["track"].unique()):
        dtrack = df[df["track"] == track]
        for N in sorted(dtrack["N"].unique()):
            s = summarize(dtrack[dtrack["N"] == N], args.ci)
            s.to_csv(outdir / f"table_{track}_N{N}.csv", index=False)
            (outdir / f"table_{track}_N{N}.tex").write_text(
                to_latex(s, f"Tree generation on held-out roots ({track}, N={N})"))
        main = summarize(dtrack, args.ci)
        main.to_csv(outdir / f"table_{track}_main.csv", index=False)
        (outdir / f"table_{track}_main.tex").write_text(
            to_latex(main, f"Tree generation on held-out roots ({track}, all N)"))
        print(f"[{track}] wrote tables for N in {sorted(dtrack['N'].unique())} + main")


if __name__ == "__main__":
    main()
