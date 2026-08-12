"""Metric C — amino-acid substitution spectrum."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from antibody_benchmark.metrics._stats import bootstrap_ci, js_divergence, spearman
from antibody_benchmark.trees import AntibodyTree

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}


def _edge_subs(tree: AntibodyTree, seqs: Mapping[str, str]) -> np.ndarray:
    M = np.zeros((20, 20), dtype=float)
    for p, c in tree.edges:
        sp, sc = seqs[p], seqs[c]
        if len(sp) != len(sc):
            continue
        for a, b in zip(sp, sc):
            if a == b:
                continue
            if a in AA_IDX and b in AA_IDX:
                M[AA_IDX[a], AA_IDX[b]] += 1.0
    return M


def metric_substitution_spectrum(
    trees: Sequence[AntibodyTree],
    generated_by_family: Mapping[str, Sequence[Mapping[str, str]]],
    *,
    n_boot: int = 1000,
    seed: int = 0,
) -> dict:
    jss, rhos = [], []
    for tree in trees:
        real_seqs = {
            n: tree.true_aa_sequences.get(n) or tree.true_sequences[n] for n in tree.node_ids
        }
        Mr = _edge_subs(tree, real_seqs)
        gens = generated_by_family.get(tree.family_id, [])
        if not gens or Mr.sum() == 0:
            continue
        Mg = sum((_edge_subs(tree, g) for g in gens), start=np.zeros_like(Mr))
        # off-diagonal vector
        mask = ~np.eye(20, dtype=bool)
        rv, gv = Mr[mask], Mg[mask]
        jss.append(js_divergence(rv, gv))
        rhos.append(spearman(rv, gv))
    return {
        "metric": "C_substitution_spectrum",
        "js_divergence": bootstrap_ci(jss, n_boot=n_boot, seed=seed),
        "offdiag_spearman": bootstrap_ci(rhos, n_boot=n_boot, seed=seed + 1),
        "n_families": len(jss),
        "status": "ok" if jss else "empty",
    }
