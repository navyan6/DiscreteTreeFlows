"""
Fit birth–death + substitution parameters on TRAIN trees only (no leakage).

Reconstructed trees carry no extinct lineages, so death is not identifiable —
we report a Yule (pure-birth) fit (death=0) and say so. Birth rate via a moment
estimator (b ≈ (n_tips-1)/total_tree_length), amino-acid frequencies empirical.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

from benchmarks.heldout.build_examples import load_tree, list_groups
from benchmarks.metrics import trees as T

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"


def fit_params(train_dir: str | Path) -> dict:
    births, all_bl = [], []
    freqs = Counter()
    for g in list_groups(train_dir):
        t = load_tree(train_dir, g)
        leaves = T.leaf_labels(t)
        total_len = sum(t.branch_lengths.values())
        if total_len > 0 and len(leaves) > 1:
            births.append((len(leaves) - 1) / total_len)
        all_bl.extend(t.branch_lengths.values())
        for l in leaves:
            freqs.update(c for c in t.node_seqs.get(l, "") if c in AA_VOCAB)
    tot = sum(freqs[a] for a in AA_VOCAB) or 1
    return {
        "birth": 1.0,   # Yule shape is rate-invariant; timescale is set per-example by H
        "death": 0.0,   # death not identifiable from reconstructed (extant-only) trees
        "empirical_birth_rate": float(median(births)) if births else 1.0,  # for the record
        "subst_scale": 1.0,   # branch lengths already in subs/site
        "mean_branch_length": float(sum(all_bl) / len(all_bl)) if all_bl else 0.0,
        "aa_freqs": {a: freqs[a] / tot for a in AA_VOCAB},
        "note": "Yule pure-birth; birth normalized to 1 (topology rate-invariant, times set by H).",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-data", default="data/h3n2/train")
    ap.add_argument("--out", default="benchmarks/results/params.json")
    args = ap.parse_args()
    p = fit_params(ROOT / args.train_data)
    Path(ROOT / args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(ROOT / args.out).write_text(json.dumps(p, indent=2))
    print(json.dumps({k: v for k, v in p.items() if k != "aa_freqs"}, indent=2))


if __name__ == "__main__":
    main()
