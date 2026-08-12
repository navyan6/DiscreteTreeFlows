#!/usr/bin/env python3
"""Run root-conditioned recursive rollouts for enabled models."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antibody_benchmark.models.base import load_models
from antibody_benchmark.rollout.sample_all import sample_all
from antibody_benchmark.trees import AntibodyTree


def _load_trees(path: Path) -> list[AntibodyTree]:
    """Load frozen trees from JSONL or legacy JSON list."""
    text = path.read_text().strip()
    if not text:
        return []
    if path.suffix == ".jsonl" or text.startswith("{"):
        raw = [json.loads(line) for line in text.splitlines() if line.strip()]
    else:
        raw = json.loads(text)
    trees = []
    for d in raw:
        bl = {}
        for k, v in d["branch_lengths"].items():
            a, b = k.split("->")
            bl[(a, b)] = float(v)
        edges = [tuple(e) for e in d["edges"]]
        trees.append(
            AntibodyTree(
                family_id=d["family_id"],
                root_id=d["root_id"],
                node_ids=d.get("node_ids")
                or sorted({n for e in edges for n in e}),
                edges=edges,
                branch_lengths=bl,
                true_sequences=d.get("true_sequences")
                or {d["root_id"]: d["root_sequence"]},
                is_leaf=d.get("is_leaf")
                or {n: n in d.get("leaf_ids", []) for n in d.get("node_ids", [])},
                donor_id=d.get("donor_id", ""),
                true_aa_sequences=d.get("true_aa_sequences")
                or (
                    {d["root_id"]: d["root_aa_sequence"]}
                    if d.get("root_aa_sequence")
                    else {}
                ),
                metadata=d.get("metadata", {}),
            )
        )
    return trees


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="antibody_benchmark/configs/smoke.yaml")
    ap.add_argument(
        "--trees",
        default=None,
        help="Override trees path (default: processed/benchmark_trees.jsonl)",
    )
    ap.add_argument("--max-families", type=int, default=None)
    ap.add_argument("--n-rollouts", type=int, default=None)
    ap.add_argument(
        "--models",
        default=None,
        help="Comma-separated model keys to enable (overrides YAML enabled flags)",
    )
    args = ap.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text())
    if args.models:
        wanted = {m.strip() for m in args.models.split(",") if m.strip()}
        mcfg = cfg.setdefault("models", {})
        for name, entry in list(mcfg.items()):
            if not isinstance(entry, dict):
                continue
            entry["enabled"] = name in wanted
    proc = ROOT / cfg["paths"]["processed"]
    trees_path = Path(args.trees) if args.trees else proc / "benchmark_trees.jsonl"
    if not trees_path.is_file():
        legacy = proc / "heldout_trees.json"
        if legacy.is_file():
            trees_path = legacy
        else:
            raise SystemExit(f"Run build_dataset.py first; missing {trees_path}")
    trees = _load_trees(trees_path)
    if args.max_families is not None:
        trees = trees[: args.max_families]
    try:
        models = load_models(cfg)
    except RuntimeError as e:
        print(str(e), file=sys.stderr)
        raise SystemExit(2) from e
    if not models:
        raise SystemExit("No models enabled in config")

    out = ROOT / cfg["paths"]["results"] / "samples"
    n_rollouts = int(args.n_rollouts if args.n_rollouts is not None else cfg.get("n_rollouts", 20))
    records = sample_all(
        trees,
        models,
        out_dir=out,
        n_rollouts=n_rollouts,
        base_seed=int(cfg.get("seed", 42)),
        use_aa_for=frozenset({"cosine", "treesbm", "identity_null", "poisson_null"}),
    )
    print(
        json.dumps(
            {
                "n_sample_files": len(records),
                "n_families": len(trees),
                "n_rollouts": n_rollouts,
                "models": list(models),
                "trees_path": str(trees_path),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
