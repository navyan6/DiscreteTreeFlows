#!/usr/bin/env python3
"""
Export a pool of TRAIN-set topologies for ARTreeFormer/PhyloVAE training
(benchmarks/EXTERNAL.md step 1: "training topologies").

Reuses build_examples() -- same (tree, internal-root) selection and induced-
subtree sampling the held-out-root task uses for its target subtrees -- just
applied to TRAIN (not TEST) trees, with branch lengths and sequences stripped.

Leaves are anonymized to a fixed generic label set "0".."N-1" (not the real
H3N2 strain ids). Both external repos hard-require a single shared taxon set
across every tree in a training dataset (they build per-taxon identity
embeddings from one `taxa` list -- see ARTreeFormer/PhyloVAE's own
process_data()); our pool instead has many *different* subtrees with
different real leaf sets. Anonymizing makes every example share the exact
same taxon alphabet, so what gets learned is the exchangeable distribution
over N-leaf tree shapes (an unconditional topology prior), not per-strain
identity -- consistent with how TopologyPriorMethod only ever consumes bare
topology shape (leaf identity is thrown away and sequence-matched downstream
regardless, see benchmarks/methods/topology_prior.py).

Output per N:
  - train_topologies_N{N}.nwk      bare newick, one per line (reference/debug)
  - train_topologies_N{N}.trprobs  NEXUS trees block, uniform-weighted, in the
                                    exact format Bio.Phylo.parse(..., 'nexus')
                                    / mcmc_treeprob() expects (validated locally
                                    via round-trip through Phylo.write ->
                                    ete3.Tree, the same path both repos use)

Usage:
    python benchmarks/heldout/export_train_topologies.py \
        --data-dir data/h3n2/train --N 16 32 64 --per-tree 20 \
        --out-dir benchmarks/external_pools
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from benchmarks.heldout.build_examples import build_examples


def topology_newick(edges: list[tuple[str, str]], root: str, relabel: dict[str, str]) -> str:
    """Bare topology newick: nested parens + anonymized leaf names, no branch lengths/labels."""
    cm: dict[str, list[str]] = {}
    for p, c in edges:
        cm.setdefault(p, []).append(c)

    def render(n: str) -> str:
        kids = cm.get(n, [])
        if not kids:
            return relabel[n]
        return "(" + ",".join(render(c) for c in kids) + ")"

    return render(root) + ";"


def anonymize_leaves(edges: list[tuple[str, str]], root: str, n: int) -> dict[str, str]:
    """Map this example's actual leaf ids -> the fixed generic alphabet "0".."n-1",
    in the order leaves appear via preorder traversal (order is arbitrary but
    deterministic; what matters is every example uses the same n-leaf alphabet)."""
    cm: dict[str, list[str]] = {}
    for p, c in edges:
        cm.setdefault(p, []).append(c)
    leaves = []
    stack = [root]
    while stack:
        node = stack.pop()
        kids = cm.get(node, [])
        if not kids:
            leaves.append(node)
        else:
            stack.extend(reversed(kids))
    assert len(leaves) == n, f"expected {n} leaves, found {len(leaves)}"
    return {leaf: str(i) for i, leaf in enumerate(leaves)}


def write_trprobs(topologies: list[str], out_path: Path):
    """Uniform-weighted NEXUS trees block; literal leaf names (no Translate
    block needed -- validated against Bio.Phylo.parse(..., 'nexus'))."""
    n = len(topologies)
    w = 1.0 / n
    lines = ["#NEXUS", "Begin trees;"]
    for i, nwk in enumerate(topologies, start=1):
        lines.append(f"    tree tree_{i} = [&W {w:.10f}] [&U] {nwk}")
    lines.append("End;")
    out_path.write_text("\n".join(lines) + "\n")


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
            relabel = anonymize_leaves(edges, ex["root_id"], N)
            lines.append(topology_newick(edges, ex["root_id"], relabel))

        nwk_path = out_dir / f"train_topologies_N{N}.nwk"
        nwk_path.write_text("\n".join(lines) + "\n")

        trprobs_path = out_dir / f"train_topologies_N{N}.trprobs"
        write_trprobs(lines, trprobs_path)

        print(f"N={N}: {len(lines)} topologies from {len({e['group'] for e in examples})} "
              f"train trees -> {nwk_path.name}, {trprobs_path.name}")


if __name__ == "__main__":
    main()
