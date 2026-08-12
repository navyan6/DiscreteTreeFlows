"""Metric A — root-to-leaf mutation-distance distribution."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from antibody_benchmark.metrics._stats import bootstrap_ci, spearman, wasserstein1
from antibody_benchmark.trees import AntibodyTree, hamming


def _aa(tree: AntibodyTree, node: str, seqs: Mapping[str, str] | None = None) -> str:
    if seqs is not None:
        return seqs[node]
    return tree.true_aa_sequences.get(node) or tree.true_sequences[node]


def root_to_leaf_distances(tree: AntibodyTree, seqs: Mapping[str, str] | None = None) -> list[float]:
    root = _aa(tree, tree.root_id, seqs)
    out = []
    for leaf in tree.leaves():
        leaf_seq = _aa(tree, leaf, seqs)
        if len(root) == len(leaf_seq):
            out.append(float(hamming(root, leaf_seq)))
    return out


def metric_root_to_leaf(
    trees: Sequence[AntibodyTree],
    generated_by_family: Mapping[str, Sequence[Mapping[str, str]]],
    *,
    n_boot: int = 1000,
    seed: int = 0,
) -> dict:
    """Per-family W1 + mean-error; aggregate with family bootstrap."""
    w1s, mean_errs, real_means, gen_means = [], [], [], []
    for tree in trees:
        real = root_to_leaf_distances(tree)
        gens = generated_by_family.get(tree.family_id, [])
        if not real or not gens:
            continue
        gen_all = []
        for g in gens:
            gen_all.extend(root_to_leaf_distances(tree, g))
        w1s.append(wasserstein1(real, gen_all))
        rm, gm = float(np.mean(real)), float(np.mean(gen_all))
        real_means.append(rm)
        gen_means.append(gm)
        mean_errs.append(abs(gm - rm))
    return {
        "metric": "A_root_to_leaf",
        "wasserstein1": bootstrap_ci(w1s, n_boot=n_boot, seed=seed),
        "mean_distance_error": bootstrap_ci(mean_errs, n_boot=n_boot, seed=seed + 1),
        "mean_distance_correlation": spearman(real_means, gen_means),
        "n_families": len(w1s),
        "status": "ok" if w1s else "empty",
    }
