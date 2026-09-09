#!/usr/bin/env python3
"""
Compare VaxSeer forecasted dominance of generated flu leaves vs observed season HA.

Requires cloned wxsh1213/vaxseer + Dropbox lm weights (see adapter docstring).

Example:
  export VAXSEER_ROOT=$LABHOME/third_party/vaxseer
  python scripts/eval_vaxseer_dominance.py \
    --subtype a_h3n2 --year 2024 \
    --observed path/to/season_Y_ha.fasta \
    --generated treesbm:path/to/gen.fasta neutral:path/to/neu.fasta \
    --out benchmarks/results/tables/vaxseer_dominance_h3n2_2024.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from benchmarks.methods.vaxseer_dominance import (
    VaxSeerDominanceScorer,
    compare_generated_vs_observed,
)


def read_fasta(path: Path) -> list[str]:
    seqs, cur = [], []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper().replace("-", ""))
                cur = []
            else:
                cur.append(line)
        if cur:
            seqs.append("".join(cur).upper().replace("-", ""))
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subtype", default="a_h3n2", choices=["a_h3n2", "a_h1n1"])
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--observed", type=Path, required=True)
    ap.add_argument(
        "--generated",
        nargs="+",
        required=True,
        help="method:fasta pairs, e.g. treesbm:gen.fa neutral:n.fa",
    )
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--md-out", type=Path, default=None)
    args = ap.parse_args()

    observed = read_fasta(args.observed)
    generated = {}
    for item in args.generated:
        if ":" not in item:
            raise SystemExit(f"bad --generated item {item!r}; want method:path")
        method, path = item.split(":", 1)
        generated[method] = read_fasta(Path(path))

    scorer = VaxSeerDominanceScorer(
        subtype=args.subtype, year=args.year, device=args.device
    )
    result = compare_generated_vs_observed(generated, observed, scorer=scorer)
    result["subtype"] = args.subtype
    result["year"] = args.year
    result["observed_path"] = str(args.observed)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {args.out}")

    md = args.md_out or args.out.with_suffix(".md")
    lines = [
        f"# VaxSeer dominance — {args.subtype} {args.year}",
        "",
        f"Observed n={result['n_observed']}  top-decile n={result['n_top_decile']}",
        f"p_obs mean/median/max = {result['p_obs_mean']:.4g} / "
        f"{result['p_obs_median']:.4g} / {result['p_obs_max']:.4g}",
        "",
        "| method | n | p_mean | p_median | p_max | top10%@ε2 | top10%@ε5 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for m, d in result["methods"].items():
        lines.append(
            f"| {m} | {d['n_generated']} | {d['p_mean']:.4g} | {d['p_median']:.4g} | "
            f"{d['p_max']:.4g} | {d.get('top_decile_recovery_eps2')} | "
            f"{d.get('top_decile_recovery_eps5')} |"
        )
    lines += [
        "",
        "TreeSBM = generator; VaxSeer lm = frozen external scorer.",
        "See benchmarks/results/tables/vaxseer_dominance.md.",
        "",
    ]
    md.write_text("\n".join(lines) + "\n")
    print(f"Wrote {md}")


if __name__ == "__main__":
    main()
