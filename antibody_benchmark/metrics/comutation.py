"""Metric D — pairwise co-mutation / epistatic structure."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from antibody_benchmark.metrics._stats import bootstrap_ci, spearman
from antibody_benchmark.trees import AntibodyTree


def _comutation_matrix(tree: AntibodyTree, seqs_list: Sequence[Mapping[str, str]], min_freq: float = 0.05) -> np.ndarray:
    root_id = tree.root_id
    L = None
    mut_counts = None
    pair_counts = None
    n = 0
    for seqs in seqs_list:
        root = seqs[root_id]
        if L is None:
            L = len(root)
            mut_counts = np.zeros(L)
            pair_counts = np.zeros((L, L))
        for leaf in tree.leaves():
            s = seqs[leaf]
            if len(s) != L:
                continue
            mut = np.fromiter((a != b for a, b in zip(root, s)), dtype=float, count=L)
            mut_counts += mut
            pair_counts += np.outer(mut, mut)
            n += 1
    if not n or mut_counts is None:
        return np.zeros((0, 0))
    p = mut_counts / n
    C = pair_counts / n
    # C_ij / (p_i p_j); diagonal ignored later
    denom = np.outer(p, p)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(denom > 0, C / denom, 0.0)
    # restrict to frequent sites
    keep = p >= min_freq
    return out[np.ix_(keep, keep)]


def metric_comutation(
    trees: Sequence[AntibodyTree],
    generated_by_family: Mapping[str, Sequence[Mapping[str, str]]],
    *,
    n_boot: int = 1000,
    seed: int = 0,
) -> dict:
    rhos, frobs = [], []
    for tree in trees:
        real = _comutation_matrix(
            tree,
            [{n: tree.true_aa_sequences.get(n) or tree.true_sequences[n] for n in tree.node_ids}],
        )
        gens = generated_by_family.get(tree.family_id, [])
        if real.size == 0 or not gens:
            continue
        gen = _comutation_matrix(tree, gens)
        m = min(real.shape[0], gen.shape[0])
        if m < 2:
            continue
        iu = np.triu_indices(m, k=1)
        rv, gv = real[:m, :m][iu], gen[:m, :m][iu]
        rhos.append(spearman(rv, gv))
        denom = np.linalg.norm(rv) + 1e-12
        frobs.append(float(np.linalg.norm(rv - gv) / denom))
    return {
        "metric": "D_comutation",
        "spearman": bootstrap_ci(rhos, n_boot=n_boot, seed=seed),
        "normalized_frobenius": bootstrap_ci(frobs, n_boot=n_boot, seed=seed + 1),
        "n_families": len(rhos),
        "status": "ok" if rhos else "empty",
    }
