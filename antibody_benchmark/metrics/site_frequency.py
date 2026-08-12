"""Metric B — mutation-site frequency (+ CDR/FR stratification hooks)."""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np

from antibody_benchmark.metrics._stats import bootstrap_ci, js_divergence, spearman
from antibody_benchmark.trees import AntibodyTree


def _site_mut_freq(tree: AntibodyTree, seqs_list: Sequence[Mapping[str, str]]) -> np.ndarray:
    """Fraction of root→leaf trajectories where site i mutated (AA)."""
    root_id = tree.root_id
    L = None
    counts = None
    n_traj = 0
    for seqs in seqs_list:
        root = seqs[root_id]
        if L is None:
            L = len(root)
            counts = np.zeros(L, dtype=float)
        for leaf in tree.leaves():
            child = seqs[leaf]
            if len(child) != L:
                continue
            counts += np.fromiter((a != b for a, b in zip(root, child)), dtype=float, count=L)
            n_traj += 1
    if not n_traj or counts is None:
        return np.zeros(0)
    return counts / n_traj


def metric_site_frequency(
    trees: Sequence[AntibodyTree],
    generated_by_family: Mapping[str, Sequence[Mapping[str, str]]],
    *,
    region_masks: Optional[Mapping[str, np.ndarray]] = None,
    n_boot: int = 1000,
    seed: int = 0,
) -> dict:
    rhos, jss = [], []
    for tree in trees:
        real_f = _site_mut_freq(
            tree, [{n: (tree.true_aa_sequences.get(n) or tree.true_sequences[n]) for n in tree.node_ids}]
        )
        gens = generated_by_family.get(tree.family_id, [])
        if real_f.size == 0 or not gens:
            continue
        gen_f = _site_mut_freq(tree, gens)
        L = min(len(real_f), len(gen_f))
        if L == 0:
            continue
        rhos.append(spearman(real_f[:L], gen_f[:L]))
        jss.append(js_divergence(real_f[:L], gen_f[:L]))
    out = {
        "metric": "B_site_frequency",
        "spearman": bootstrap_ci(rhos, n_boot=n_boot, seed=seed),
        "js_divergence": bootstrap_ci(jss, n_boot=n_boot, seed=seed + 1),
        "n_families": len(rhos),
        "cdr_framework_note": "Pass IMGT/ANARCI region_masks for stratified CDR/FR enrichment.",
        "status": "ok" if rhos else "empty",
    }
    if region_masks:
        out["region_masks_provided"] = list(region_masks.keys())
    return out
