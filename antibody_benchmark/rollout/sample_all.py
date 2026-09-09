"""Run N rollouts for all trees × models and persist samples."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Mapping

from antibody_benchmark.models.base import EvolutionModel
from antibody_benchmark.rollout.rollout_tree import rollout_tree
from antibody_benchmark.trees import AntibodyTree


def _root_sequence_for_model(tree: AntibodyTree, use_aa: bool) -> str:
    if use_aa:
        return tree.true_aa_sequences.get(tree.root_id) or tree.true_sequences[tree.root_id]
    return tree.true_sequences[tree.root_id]


def sample_all(
    trees: Iterable[AntibodyTree],
    models: Mapping[str, EvolutionModel],
    *,
    out_dir: Path,
    n_rollouts: int = 20,
    base_seed: int = 42,
    use_aa_for: frozenset[str] | None = None,
) -> list[dict]:
    """Write JSON samples under out_dir/<model>/<family>/rollout_k.json.

    Rollout receives only root_sequence + edges + branch_lengths (never
    non-root true_sequences).
    """
    use_aa_for = use_aa_for or frozenset({"cosine", "treesbm", "identity_null"})
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []

    for model_name, model in models.items():
        for tree in trees:
            fam_dir = out_dir / model_name / tree.family_id.replace("/", "_")
            fam_dir.mkdir(parents=True, exist_ok=True)
            use_aa = model_name in use_aa_for or getattr(model, "alphabet", "nt") == "aa"
            root_seq = _root_sequence_for_model(tree, use_aa)
            for k in range(n_rollouts):
                seqs = rollout_tree(
                    model,
                    root_seq,
                    tree.edges,
                    tree.branch_lengths,
                    seed=base_seed + 1_000_003 * k,
                    root_id=tree.root_id,
                    family_id=tree.family_id,
                )
                path = fam_dir / f"rollout_{k:03d}.json"
                payload = {
                    "family_id": tree.family_id,
                    "model": model_name,
                    "rollout_idx": k,
                    "root_id": tree.root_id,
                    "sequences": seqs,
                    "alphabet": "aa" if use_aa else "nt",
                }
                path.write_text(json.dumps(payload))
                records.append(
                    {
                        "family_id": tree.family_id,
                        "model": model_name,
                        "rollout_idx": k,
                        "path": str(path),
                    }
                )
    manifest = out_dir / "manifest.json"
    manifest.write_text(json.dumps(records, indent=2))
    return records
