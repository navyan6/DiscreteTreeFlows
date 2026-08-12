#!/usr/bin/env python3
"""Evaluate saved rollouts; write primary/secondary TEMP tables."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antibody_benchmark.codon import looks_like_nt, translate_nt
from antibody_benchmark.metrics.comutation import metric_comutation
from antibody_benchmark.metrics.diversity import metric_leaf_diversity
from antibody_benchmark.metrics.likelihood import exact_descendant_hamming
from antibody_benchmark.metrics.mutation_distance import metric_root_to_leaf
from antibody_benchmark.metrics.site_frequency import metric_site_frequency
from antibody_benchmark.metrics.substitution_spectrum import metric_substitution_spectrum
from antibody_benchmark.metrics.tree_metrics import metric_tree_shape_scaffold
from antibody_benchmark.trees import AntibodyTree


def _load_trees(path: Path) -> list[AntibodyTree]:
    text = path.read_text().strip()
    if path.suffix == ".jsonl" or (text and text.startswith("{")):
        raw = [json.loads(line) for line in text.splitlines() if line.strip()]
    else:
        raw = json.loads(text)
    trees = []
    for d in raw:
        bl = {tuple(k.split("->")): float(v) for k, v in d["branch_lengths"].items()}
        edges = [tuple(e) for e in d["edges"]]
        trees.append(
            AntibodyTree(
                family_id=d["family_id"],
                root_id=d["root_id"],
                node_ids=d.get("node_ids") or sorted({n for e in edges for n in e}),
                edges=edges,
                branch_lengths=bl,
                true_sequences=d.get("true_sequences")
                or {d["root_id"]: d["root_sequence"]},
                is_leaf=d.get("is_leaf")
                or {n: n in set(d.get("leaf_ids", [])) for n in d.get("node_ids", [])},
                donor_id=d.get("donor_id", ""),
                true_aa_sequences=d.get("true_aa_sequences", {}),
                metadata=d.get("metadata", {}),
            )
        )
    return trees


def _to_aa(seqs: dict[str, str]) -> dict[str, str]:
    out = {}
    for k, v in seqs.items():
        out[k] = translate_nt(v) if looks_like_nt(v) else v
    return out


def _load_generated(samples_dir: Path, model: str) -> dict[str, list[dict[str, str]]]:
    by_fam: dict[str, list[dict[str, str]]] = defaultdict(list)
    model_dir = samples_dir / model
    if not model_dir.is_dir():
        return {}
    for path in sorted(model_dir.glob("*/*.json")):
        payload = json.loads(path.read_text())
        seqs = _to_aa(payload["sequences"])
        by_fam[payload["family_id"]].append(seqs)
    return by_fam


def _fmt_ci(d: dict | None) -> str:
    if not d or d.get("point") != d.get("point"):
        return "TEMP/—"
    return f"{d['point']:.3f} [{d['lo']:.3f},{d['hi']:.3f}]"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="antibody_benchmark/configs/smoke.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text())
    proc = ROOT / cfg["paths"]["processed"]
    trees_path = proc / "benchmark_trees.jsonl"
    if not trees_path.is_file():
        trees_path = proc / "heldout_trees.json"
    trees = _load_trees(trees_path)
    samples_dir = ROOT / cfg["paths"]["results"] / "samples"
    summary_dir = ROOT / cfg["paths"]["results"] / "summary"
    per_tree_dir = ROOT / cfg["paths"]["results"] / "per_tree"
    summary_dir.mkdir(parents=True, exist_ok=True)
    per_tree_dir.mkdir(parents=True, exist_ok=True)

    model_names = sorted({p.parent.parent.name for p in samples_dir.glob("*/*/rollout_*.json")})
    # Always include paper rows even if empty
    paper_models = ["thrifty", "dasm_thrifty", "cosine", "treesbm"]
    for m in paper_models:
        if m not in model_names:
            model_names.append(m)

    primary_rows = []
    secondary_rows = []
    detailed = {}
    n_boot = int(cfg.get("n_bootstrap", 1000))
    seed = int(cfg.get("seed", 0))

    for model in model_names:
        gens = _load_generated(samples_dir, model)
        stub = not gens
        label_suffix = " TEMP/stub" if stub or model in {"identity_null", "poisson_null"} else ""
        if stub:
            primary_rows.append(
                {
                    "Model": model + label_suffix,
                    "Root→leaf W1 ↓": "TEMP/—",
                    "Site freq ρ ↑": "TEMP/—",
                    "Substitution JS ↓": "TEMP/—",
                    "Co-mutation ρ ↑": "TEMP/—",
                    "Leaf diversity W1 ↓": "TEMP/—",
                }
            )
            secondary_rows.append(
                {
                    "Model": model + label_suffix,
                    "Transition NLL ↓": "N/A" if model == "treesbm" else "TEMP/—",
                    "Branch mutation calibration ↓": "TEMP/—",
                    "Exact descendant distance ↓": "TEMP/—",
                }
            )
            continue

        A = metric_root_to_leaf(trees, gens, n_boot=n_boot, seed=seed)
        B = metric_site_frequency(trees, gens, n_boot=n_boot, seed=seed)
        C = metric_substitution_spectrum(trees, gens, n_boot=n_boot, seed=seed)
        D = metric_comutation(trees, gens, n_boot=n_boot, seed=seed)
        E = metric_leaf_diversity(trees, gens, n_boot=n_boot, seed=seed)
        exact = exact_descendant_hamming(trees, gens)
        detailed[model] = {"A": A, "B": B, "C": C, "D": D, "E": E, "exact": exact}
        primary_rows.append(
            {
                "Model": model + label_suffix,
                "Root→leaf W1 ↓": _fmt_ci(A.get("wasserstein1")),
                "Site freq ρ ↑": _fmt_ci(B.get("spearman")),
                "Substitution JS ↓": _fmt_ci(C.get("js_divergence")),
                "Co-mutation ρ ↑": _fmt_ci(D.get("spearman")),
                "Leaf diversity W1 ↓": _fmt_ci(E.get("wasserstein1")),
            }
        )
        secondary_rows.append(
            {
                "Model": model + label_suffix,
                "Transition NLL ↓": "N/A" if model == "treesbm" else "TEMP/—",
                "Branch mutation calibration ↓": "TEMP/—",
                "Exact descendant distance ↓": (
                    f"{exact['mean_frac_hamming']:.3f}"
                    if exact.get("mean_frac_hamming") is not None
                    else "TEMP/—"
                ),
            }
        )

    primary = pd.DataFrame(primary_rows)
    secondary = pd.DataFrame(secondary_rows)
    primary.to_csv(summary_dir / "primary_metrics.csv", index=False)
    secondary.to_csv(summary_dir / "secondary_metrics.csv", index=False)
    (summary_dir / "detailed_metrics.json").write_text(json.dumps(detailed, indent=2, default=str))
    (summary_dir / "tree_generation_scaffold.json").write_text(
        json.dumps(metric_tree_shape_scaffold(None), indent=2)
    )

    def _df_to_md(df: pd.DataFrame) -> str:
        cols = list(df.columns)
        lines = [
            "| " + " | ".join(cols) + " |",
            "| " + " | ".join("---" for _ in cols) + " |",
        ]
        for _, row in df.iterrows():
            lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
        return "\n".join(lines)

    # Markdown summary
    md = ["# Results summary (TEMP until full model rollouts)\n"]
    md.append("## Primary (distributional rollout)\n")
    md.append(_df_to_md(primary))
    md.append("\n\n## Secondary\n")
    md.append(_df_to_md(secondary))
    md.append(
        "\n\nNotes:\n"
        "- Recursive root-conditioned rollouts; generated parents feed children.\n"
        "- Family-level bootstrap CIs when computed.\n"
        "- Rows labeled TEMP/stub are null models or missing adapters.\n"
        "- CoSiNE must use unguided Gillespie; TreeSBM forced to observed topology+BL.\n"
    )
    (summary_dir / "RESULTS.md").write_text("\n".join(md))
    print(primary.to_string(index=False))
    print(f"\nWrote {summary_dir / 'primary_metrics.csv'}")
    print(f"Wrote {summary_dir / 'RESULTS.md'}")


if __name__ == "__main__":
    main()
