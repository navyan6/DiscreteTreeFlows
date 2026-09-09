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


def _mean_root_tip_height(tree: AntibodyTree) -> float:
    """Mean path length root→leaf from observed branch lengths."""
    from collections import defaultdict, deque

    children: dict[str, list[str]] = defaultdict(list)
    for p, c in tree.edges:
        children[p].append(c)
    dist = {tree.root_id: 0.0}
    q = deque([tree.root_id])
    while q:
        u = q.popleft()
        for v in children.get(u, []):
            dist[v] = dist[u] + float(tree.branch_lengths[(u, v)])
            q.append(v)
    tips = [
        dist[n]
        for n, is_leaf in (tree.is_leaf or {}).items()
        if is_leaf and n in dist and n != tree.root_id
    ]
    if not tips:
        tips = [v for n, v in dist.items() if n != tree.root_id and n not in children]
    return sum(tips) / len(tips) if tips else 0.05


def _n_leaves(tree: AntibodyTree) -> int:
    if tree.is_leaf:
        return sum(1 for n, v in tree.is_leaf.items() if v and n != tree.root_id)
    return max(0, len(tree.node_ids) - 1)


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
    non-root true_sequences). Free-topology models (``free_topology=True``)
    use ``sample_lineage`` instead.
    """
    use_aa_for = use_aa_for or frozenset(
        {"cosine", "treesbm", "identity_null", "plm_prior", "ar_tree_edit"}
    )
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []

    for model_name, model in models.items():
        for tree in trees:
            fam_dir = out_dir / model_name / tree.family_id.replace("/", "_")
            fam_dir.mkdir(parents=True, exist_ok=True)
            use_aa = model_name in use_aa_for or getattr(model, "alphabet", "nt") == "aa"
            root_seq = _root_sequence_for_model(tree, use_aa)
            free = bool(getattr(model, "free_topology", False))
            for k in range(n_rollouts):
                path = fam_dir / f"rollout_{k:03d}.json"
                # Resume-friendly: keep existing rollouts when extending n_rollouts.
                if path.is_file():
                    records.append(
                        {
                            "family_id": tree.family_id,
                            "model": model_name,
                            "rollout_idx": k,
                            "path": str(path),
                            "skipped_existing": True,
                        }
                    )
                    continue
                seed = base_seed + 1_000_003 * k + abs(hash(tree.family_id)) % 10_007
                if free and hasattr(model, "sample_lineage"):
                    seqs = model.sample_lineage(  # type: ignore[attr-defined]
                        root_seq,
                        _n_leaves(tree),
                        _mean_root_tip_height(tree),
                        seed,
                    )
                    root_id_out = tree.root_id
                    for nid, s in seqs.items():
                        if s == root_seq:
                            root_id_out = nid
                            break
                else:
                    seqs = rollout_tree(
                        model,
                        root_seq,
                        tree.edges,
                        tree.branch_lengths,
                        seed=seed,
                        root_id=tree.root_id,
                        family_id=tree.family_id,
                    )
                    root_id_out = tree.root_id
                payload = {
                    "family_id": tree.family_id,
                    "model": model_name,
                    "rollout_idx": k,
                    "root_id": root_id_out,
                    "sequences": seqs,
                    "alphabet": "aa" if use_aa else "nt",
                    "free_topology": free,
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
