#!/usr/bin/env python3
"""
Export a pool of TRAIN-set topologies for ARTreeFormer/PhyloVAE training
(benchmarks/EXTERNAL.md step 1: "training topologies").

Reuses build_examples() -- same (tree, internal-root) selection and induced-
subtree sampling the held-out-root task uses for its target subtrees -- just
applied to TRAIN (not TEST) trees, with branch lengths and sequences stripped.
These external repos train tree-topology density models on bare newick sets
(bipartition structure only), not branch lengths or sequences.

Usage:
    python benchmarks/heldout/export_train_topologies.py \
        --data-dir data/h3n2/train --N 16 32 64 --per-tree 20 \
        --out-dir benchmarks/external_pools
    # -> benchmarks/external_pools/train_topologies_N16.nwk (+ N32, N64)
    # one bare newick per line, e.g. (l1,l2,(l3,l4));
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from benchmarks.heldout.build_examples import build_examples


def topology_newick(edges: list[tuple[str, str]], root: str) -> str:
    """Bare topology newick: nested parens + leaf names, no branch lengths/labels."""
    cm: dict[str, list[str]] = {}
    for p, c in edges:
        cm.setdefault(p, []).append(c)

    def render(n: str) -> str:
        kids = cm.get(n, [])
        if not kids:
            return n
        return "(" + ",".join(render(c) for c in kids) + ")"

    return render(root) + ";"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/h3n2/train")
    ap.add_argument("--N", type=int, nargs="+", default=[16, 32, 64])
    ap.add_argument("--per-tree", type=int, default=20,
                    help="Max sampled roots per tree per N (build_examples' max_roots_per_tree)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-dir", default="benchmarks/external_pools")
    args = ap.parse_args()

    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    for N in args.N:
        examples = build_examples(ROOT / args.data_dir, N, seed=args.seed,
                                  max_roots_per_tree=args.per_tree)
        lines = []
        for ex in examples:
            edges = [tuple(k.split("|")) for k in ex["target_branch_lengths"]]
            lines.append(topology_newick(edges, ex["root_id"]))
        out_path = out_dir / f"train_topologies_N{N}.nwk"
        out_path.write_text("\n".join(lines) + "\n")
        print(f"N={N}: {len(lines)} topologies from {len({e['group'] for e in examples})} "
              f"train trees -> {out_path}")


if __name__ == "__main__":
    main()
