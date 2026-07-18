"""
Held-out-root example construction (Task definition + Dataset construction).

For each processed tree we pick internal nodes as roots, extract the descendant
subtree, select exactly N terminal descendants, take the induced subtree,
collapse degree-2 nodes (summing branch lengths), and record the horizon H (mean
root-to-tip patristic distance). Roots come from TEST trees; H buckets and all
fitted parameters use TRAIN trees only. Train/test are disjoint trees (temporal
split), so there is no node leakage; we assert it.

Each example gives a method only {root_seq, N, H}; the induced subtree (topology
+ branch lengths + terminal/internal sequences) is the hidden target.

Pure-Python (Bio.Phylo + the metric library); no torch/model dependency.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import deque
from pathlib import Path

from Bio import Phylo, SeqIO

ROOT = Path(__file__).resolve().parents[2]
import sys
sys.path.insert(0, str(ROOT))

from src.tree_state import TreeState
from benchmarks.metrics import trees as T


# ── loading a processed tree into a TreeState ────────────────────────────────

def _parse_newick(nwk_path: str):
    """(root_id, edges, branch_lengths) with internal nodes named NODE_*."""
    tree = Phylo.read(nwk_path, "newick")
    edges, bls, counter = [], {}, [0]

    def name(clade):
        if clade.name:
            return clade.name
        clade.name = f"NODE_{counter[0]:07d}"
        counter[0] += 1
        return clade.name

    def walk(parent):
        pn = name(parent)
        for child in parent.clades:
            cn = name(child)
            edges.append((pn, cn))
            bls[(pn, cn)] = child.branch_length or 0.0
            walk(child)

    walk(tree.root)
    return name(tree.root), edges, bls


def load_tree(data_dir: str | Path, g: int) -> TreeState:
    """Load group_{g}_rooted.nwk + group_{g}_anc_aa.fasta into a TreeState."""
    d = Path(data_dir)
    root_id, edges, bls = _parse_newick(str(d / f"group_{g:03d}_rooted.nwk"))
    seqs = {rec.id: str(rec.seq) for rec in
            SeqIO.parse(d / f"group_{g:03d}_anc_aa.fasta", "fasta")}
    node_ids = [root_id]
    seen = {root_id}
    cm = T.children_map(TreeState(  # cheap: build children map from edges
        node_ids=list({n for e in edges for n in e} | {root_id}),
        root_id=root_id, edges=edges, branch_lengths=bls,
        node_seqs={}, active_leaves=[]))
    q = deque([root_id])
    while q:
        n = q.popleft()
        for c in cm.get(n, []):
            if c not in seen:
                seen.add(c); node_ids.append(c); q.append(c)
    # fill missing seqs (some internal nodes may be absent) with gaps
    if seqs:
        ref_len = len(next(iter(seqs.values())))
        for n in node_ids:
            seqs.setdefault(n, "-" * ref_len)
    leaves = [n for n in node_ids if n not in cm]
    return TreeState(node_ids=node_ids, root_id=root_id, edges=edges,
                     branch_lengths=bls, node_seqs=seqs, active_leaves=leaves)


def list_groups(data_dir: str | Path) -> list[int]:
    d = Path(data_dir)
    return sorted(int(p.stem.split("_")[1]) for p in d.glob("group_*_rooted.nwk")
                  if (d / p.name.replace("_rooted.nwk", "_anc_aa.fasta")).exists())


# ── induced subtree with degree-2 collapse ───────────────────────────────────

def induced_subtree(tree: TreeState, root: str, leaves: list[str]) -> TreeState:
    """
    Minimal subtree of `root` spanning `leaves`, with degree-2 internal nodes
    collapsed and their branch lengths summed. `root` is always kept (it may end
    up with a single child). Sequences carried over from `tree`.
    """
    pm = T.parent_map(tree)
    leaf_set = set(leaves)

    # kept = union of root->leaf paths
    kept = {root}
    for lf in leaves:
        n = lf
        while n != root:
            kept.add(n)
            n = pm[n]

    # induced parent map + branch length (child -> parent, child -> bl)
    ipar: dict[str, str] = {}
    ibl: dict[str, float] = {}
    for n in kept:
        if n == root:
            continue
        ipar[n] = pm[n]
        ibl[n] = float(tree.branch_lengths.get((pm[n], n), 0.0))

    def children_of():
        cm: dict[str, list[str]] = {}
        for c, p in ipar.items():
            cm.setdefault(p, []).append(c)
        return cm

    # collapse pass-through nodes (single induced child, not root, not selected leaf)
    changed = True
    while changed:
        changed = False
        cm = children_of()
        for n, kids in cm.items():
            if n == root or n in leaf_set:
                continue
            if len(kids) == 1:
                c = kids[0]
                ibl[c] = ibl[c] + ibl[n]      # sum branch lengths
                ipar[c] = ipar[n]
                del ipar[n]
                del ibl[n]
                changed = True
                break

    # assemble TreeState (BFS order from root)
    edges = [(ipar[c], c) for c in ipar]
    branch_lengths = {(ipar[c], c): ibl[c] for c in ipar}
    cm = {}
    for p, c in edges:
        cm.setdefault(p, []).append(c)
    node_ids, seen, q = [root], {root}, deque([root])
    while q:
        n = q.popleft()
        for c in cm.get(n, []):
            if c not in seen:
                seen.add(c); node_ids.append(c); q.append(c)
    node_seqs = {n: tree.node_seqs.get(n, "") for n in node_ids}
    leaves_final = [n for n in node_ids if n not in cm]
    return TreeState(node_ids=node_ids, root_id=root, edges=edges,
                     branch_lengths=branch_lengths, node_seqs=node_seqs,
                     active_leaves=leaves_final)


def root_to_tip(tree: TreeState) -> list[float]:
    times = T.node_times(tree)
    return [times[l] for l in T.leaf_labels(tree)]


# ── example construction ─────────────────────────────────────────────────────

def build_examples(data_dir: str | Path, N: int, seed: int = 0,
                   max_roots_per_tree: int = 5, min_root_leaves: int | None = None
                   ) -> list[dict]:
    """
    One example per (tree, internal root). Returns list of dicts with root_seq,
    N, H, and the hidden target subtree (newick + node_seqs).
    """
    rng = random.Random(seed)
    min_leaves = min_root_leaves or N
    examples = []
    for g in list_groups(data_dir):
        tree = load_tree(data_dir, g)
        dl = T.descendant_leaves(tree)
        cm = T.children_map(tree)
        # candidate internal roots with enough descendant leaves
        cands = [n for n in tree.node_ids if n in cm and len(dl[n]) >= min_leaves]
        rng.shuffle(cands)
        for root in cands[:max_roots_per_tree]:
            leaves = rng.sample(sorted(dl[root]), N)
            sub = induced_subtree(tree, root, leaves)
            H = sum(root_to_tip(sub)) / N
            examples.append({
                "group": g,
                "root_id": root,
                "root_seq": tree.node_seqs[root],
                "N": N,
                "H": H,
                "target_newick": T.to_newick(sub),
                "target_seqs": {n: sub.node_seqs[n] for n in T.leaf_labels(sub)},
                "target_node_seqs": sub.node_seqs,
                "target_branch_lengths": {f"{p}|{c}": v
                                          for (p, c), v in sub.branch_lengths.items()},
            })
    return examples


def h_buckets(train_H: list[float], q=(1 / 3, 2 / 3)) -> tuple[float, float]:
    """Short/medium/long thresholds from train-set root-to-tip quantiles."""
    import numpy as np
    return tuple(float(np.quantile(train_H, x)) for x in q)


def assert_no_leakage(train_dir, test_dir):
    """
    No LEAF (real sequence id, e.g. EPI_ISL_*) may appear in both train and test.
    Internal nodes are named NODE_* fresh per tree by the parser, so they collide
    trivially and must be excluded — only leaf identities are meaningful.
    """
    def all_leaf_ids(d):
        s = set()
        for g in list_groups(d):
            root_id, edges, _ = _parse_newick(str(Path(d) / f"group_{g:03d}_rooted.nwk"))
            parents = {p for p, _ in edges}
            nodes = {n for e in edges for n in e} | {root_id}
            s |= {n for n in nodes if n not in parents}   # leaf = never a parent
        return s
    inter = all_leaf_ids(train_dir) & all_leaf_ids(test_dir)
    assert not inter, (f"LEAKAGE: {len(inter)} leaf ids shared between train and test "
                       f"(e.g. {sorted(inter)[:3]})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-data", default="data/h3n2/test")
    ap.add_argument("--train-data", default="data/h3n2/train")
    ap.add_argument("--N", type=int, nargs="+", default=[16, 32, 64])
    ap.add_argument("--max-roots-per-tree", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="benchmarks/heldout/examples")
    ap.add_argument("--check-leakage", action="store_true")
    args = ap.parse_args()

    if args.check_leakage:
        assert_no_leakage(ROOT / args.train_data, ROOT / args.test_data)
        print("leakage check passed (train ∩ test node ids = ∅)")

    outdir = ROOT / args.out
    outdir.mkdir(parents=True, exist_ok=True)
    for N in args.N:
        exs = build_examples(ROOT / args.test_data, N, seed=args.seed,
                             max_roots_per_tree=args.max_roots_per_tree)
        (outdir / f"examples_N{N}.json").write_text(json.dumps(exs))
        Hs = [e["H"] for e in exs]
        print(f"N={N}: {len(exs)} examples  H range [{min(Hs):.4f}, {max(Hs):.4f}]")


if __name__ == "__main__":
    main()
