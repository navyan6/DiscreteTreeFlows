"""Metric E — leaf repertoire diversity."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from antibody_benchmark.metrics._stats import bootstrap_ci, wasserstein1
from antibody_benchmark.trees import AntibodyTree, hamming


def _pairwise_leaf(tree: AntibodyTree, seqs: Mapping[str, str]) -> list[float]:
    leaves = tree.leaves()
    out = []
    for i in range(len(leaves)):
        for j in range(i + 1, len(leaves)):
            a, b = seqs[leaves[i]], seqs[leaves[j]]
            if len(a) == len(b):
                out.append(float(hamming(a, b)))
    return out


def metric_leaf_diversity(
    trees: Sequence[AntibodyTree],
    generated_by_family: Mapping[str, Sequence[Mapping[str, str]]],
    *,
    n_boot: int = 1000,
    seed: int = 0,
) -> dict:
    w1s, mean_errs = [], []
    for tree in trees:
        real_seqs = {
            n: tree.true_aa_sequences.get(n) or tree.true_sequences[n] for n in tree.node_ids
        }
        real = _pairwise_leaf(tree, real_seqs)
        gens = generated_by_family.get(tree.family_id, [])
        if len(real) == 0 or not gens:
            continue
        gen_all = []
        for g in gens:
            gen_all.extend(_pairwise_leaf(tree, g))
        w1s.append(wasserstein1(real, gen_all))
        mean_errs.append(abs(float(np.mean(gen_all)) - float(np.mean(real))))
    return {
        "metric": "E_leaf_diversity",
        "wasserstein1": bootstrap_ci(w1s, n_boot=n_boot, seed=seed),
        "mean_pairwise_error": bootstrap_ci(mean_errs, n_boot=n_boot, seed=seed + 1),
        "n_families": len(w1s),
        "status": "ok" if w1s else "empty",
    }
