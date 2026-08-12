#!/usr/bin/env python3
"""Incremental one-tree smoke tests for EvolutionModel adapters on Betty."""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antibody_benchmark.rollout.rollout_tree import rollout_tree
from antibody_benchmark.scripts.run_rollouts import _load_trees


def _report_path(results_dir: Path, model: str) -> Path:
    d = results_dir / "smoke_reports"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{model}.json"


def smoke_one(model_name: str, cfg: dict, tree, results_dir: Path) -> dict:
    from antibody_benchmark.models.base import load_models

    # Enable only the requested model
    mcfg = {k: {**v, "enabled": (k == model_name)} for k, v in cfg.get("models", {}).items()}
    cfg = {**cfg, "models": mcfg}
    t0 = datetime.now(timezone.utc).isoformat()
    report = {
        "model": model_name,
        "started_utc": t0,
        "family_id": tree.family_id,
        "ok": False,
        "error": None,
        "n_nodes_generated": 0,
        "root_id": tree.root_id,
    }
    try:
        models = load_models(cfg)
        if model_name not in models:
            raise RuntimeError(f"model {model_name} not loaded (enabled={mcfg.get(model_name)})")
        model = models[model_name]
        use_aa = getattr(model, "alphabet", "nt") == "aa"
        root_seq = (
            tree.true_aa_sequences.get(tree.root_id)
            if use_aa
            else tree.true_sequences[tree.root_id]
        )
        if use_aa and not root_seq:
            root_seq = tree.true_sequences[tree.root_id]
        gen = rollout_tree(
            model,
            root_seq,
            tree.edges,
            tree.branch_lengths,
            seed=0,
            root_id=tree.root_id,
            family_id=tree.family_id,
        )
        report["n_nodes_generated"] = len(gen)
        report["ok"] = tree.root_id in gen and len(gen) >= 1 + len(tree.edges)
        if not report["ok"]:
            report["error"] = f"incomplete generation: {len(gen)} nodes"
    except Exception as e:
        report["error"] = f"{type(e).__name__}: {e}"
        report["traceback"] = traceback.format_exc()[-4000:]
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    path = _report_path(results_dir, model_name)
    path.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return report


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="antibody_benchmark/configs/default.yaml")
    ap.add_argument(
        "--models",
        default="thrifty",
        help="Comma-separated: thrifty,dasm_thrifty,cosine,treesbm",
    )
    ap.add_argument("--family-index", type=int, default=0)
    args = ap.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text())
    # Betty DASM override if set in env is handled by build; trees already frozen.
    trees_path = ROOT / cfg["paths"]["processed"] / "benchmark_trees.jsonl"
    trees = _load_trees(trees_path)
    if not trees:
        raise SystemExit("No frozen trees; run build_dataset.py first")
    tree = trees[args.family_index]
    results_dir = ROOT / cfg["paths"]["results"]
    summaries = []
    for name in [m.strip() for m in args.models.split(",") if m.strip()]:
        summaries.append(smoke_one(name, cfg, tree, results_dir))
    out = results_dir / "smoke_reports" / "summary.json"
    out.write_text(json.dumps(summaries, indent=2))
    print("WROTE", out)
    if not all(s.get("ok") for s in summaries):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
