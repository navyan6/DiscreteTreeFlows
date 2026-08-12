"""Secondary: transition NLL when a model exposes normalized p(child|parent,t)."""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

from antibody_benchmark.trees import AntibodyTree


def metric_transition_nll(
    trees: Sequence[AntibodyTree],
    model: Any,
    *,
    model_name: str,
) -> dict:
    """Compute mean -1/L log p if model implements transition_nll; else N/A."""
    if not hasattr(model, "transition_nll"):
        return {
            "metric": "secondary_transition_nll",
            "model": model_name,
            "nll": None,
            "status": "N/A",
            "note": "Model does not expose a comparable normalized transition likelihood.",
        }
    vals = []
    for tree in trees:
        for p, c in tree.edges:
            parent = tree.true_sequences[p]
            child = tree.true_sequences[c]
            bl = tree.branch_lengths[(p, c)]
            try:
                vals.append(float(model.transition_nll(parent, child, bl)))
            except Exception as e:
                return {
                    "metric": "secondary_transition_nll",
                    "model": model_name,
                    "nll": None,
                    "status": "error",
                    "note": str(e),
                }
    if not vals:
        return {
            "metric": "secondary_transition_nll",
            "model": model_name,
            "nll": None,
            "status": "empty",
        }
    import statistics as stats

    return {
        "metric": "secondary_transition_nll",
        "model": model_name,
        "nll": stats.mean(vals),
        "n_edges": len(vals),
        "status": "ok",
    }


def exact_descendant_hamming(
    trees: Sequence[AntibodyTree],
    generated_by_family: Mapping[str, Sequence[Mapping[str, str]]],
) -> dict:
    """Secondary: mean Hamming(generated node, true node) — not primary."""
    from antibody_benchmark.trees import hamming
    import statistics as stats

    dists = []
    for tree in trees:
        for gens in generated_by_family.get(tree.family_id, []):
            for n in tree.node_ids:
                if n == tree.root_id:
                    continue
                true = tree.true_aa_sequences.get(n) or tree.true_sequences[n]
                pred = gens[n]
                if len(true) == len(pred):
                    dists.append(hamming(true, pred) / max(len(true), 1))
    return {
        "metric": "secondary_exact_descendant_distance",
        "mean_frac_hamming": stats.mean(dists) if dists else None,
        "n": len(dists),
        "status": "ok" if dists else "empty",
        "note": "Stochastic evolution — do not treat as primary quality metric.",
    }
