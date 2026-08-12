"""Scaffold for the SEPARATE TreeSBM tree-generation experiment.

Do not mix into the fixed-topology sequence rollout primary table.
"""

from __future__ import annotations


def metric_tree_shape_scaffold(generated_trees) -> dict:
    return {
        "metric": "tree_generation_scaffold",
        "status": "TEMP/not_run",
        "note": (
            "Evaluate TreeSBM joint generation vs BirthDeath+sequence-model separately. "
            "Metrics: normalized Sackin/Colless, cherries, depth, leaf-count, BL W1, "
            "patristic vs Hamming correlation. Do not use RF across unlabeled leaves."
        ),
        "n_trees": 0 if generated_trees is None else len(list(generated_trees)),
    }
