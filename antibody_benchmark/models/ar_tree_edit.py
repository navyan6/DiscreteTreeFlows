"""Autoregressive tree-edit baseline for Ab Track C (ARTreeFormer-adapted).

Rod.82 families have 4–12 leaves; viral ARTreeFormer pools are N∈{16,32}.
We sample an N=16 ARTreeFormer topology, prune to the family's leaf count,
rescale BLs to the observed mean root→tip height, and fill sequences with
shared JTT (pyvolve). Honest adapted-prior recipe — not Ab-trained ARTreeFormer.
"""

from __future__ import annotations

import random
from pathlib import Path

from antibody_benchmark.models.base import EvolutionModel


class ARTreeEditModel(EvolutionModel):
    """Free-topology ARTreeFormer-adapted + JTT sequences."""

    name = "ar_tree_edit"
    alphabet = "aa"
    free_topology = True

    def __init__(
        self,
        pool_path: str | None = None,
        seq_model: str = "JTT",
    ):
        repo = Path(__file__).resolve().parents[2]
        self.pool_path = Path(
            pool_path
            or repo / "benchmarks/external_pools/sampled/artreeformer_N16.nwk"
        )
        self.seq_model = seq_model
        self._pool: list[str] | None = None
        self._load_pool()

    def _load_pool(self) -> None:
        if not self.pool_path.is_file():
            raise RuntimeError(
                f"BLOCKER: missing ARTreeFormer pool at {self.pool_path}. "
                "Need benchmarks/external_pools/sampled/artreeformer_N16.nwk"
            )
        lines = [
            ln.strip()
            for ln in self.pool_path.read_text().splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        if not lines:
            raise RuntimeError(f"BLOCKER: empty ARTreeFormer pool {self.pool_path}")
        self._pool = lines

    def sample_child(self, parent_seq: str, branch_length: float, rng_seed: int) -> str:
        raise NotImplementedError(
            "ar_tree_edit uses free-topology sample_lineage(); not edge sample_child"
        )

    def sample_lineage(
        self,
        root_seq: str,
        n_leaves: int,
        tree_height: float,
        seed: int,
    ) -> dict[str, str]:
        """Return {node_id: aa_seq} for a free generated tree (root + leaves + internals)."""
        from benchmarks.adapters.sequence import evolve_pyvolve
        from benchmarks.methods.bd_topology import birth_death_topology
        from benchmarks.methods.topology_prior import parse_topology
        from benchmarks.metrics import trees as T
        from src.tree_state import TreeState

        assert self._pool is not None
        rng = random.Random(int(seed))
        n_leaves = max(int(n_leaves), 2)
        H = max(float(tree_height), 1e-6)

        topo = None
        # Prefer ARTreeFormer pool + prune to n_leaves
        for _ in range(40):
            raw = parse_topology(rng.choice(self._pool))
            leaves = T.leaf_labels(raw)
            if len(leaves) < n_leaves:
                continue
            keep = set(rng.sample(leaves, n_leaves))
            topo = _induced_subtree(raw, keep)
            if topo is not None and len(T.leaf_labels(topo)) == n_leaves:
                break
            topo = None

        if topo is None:
            # Fallback: BD topology (still generative; caption notes AR pool miss)
            topo = birth_death_topology(n_leaves, H, birth=1.0, death=0.5, seed=seed)
        else:
            # rescale BLs to mean root→tip = H
            times = T.node_times(topo)
            leaves = T.leaf_labels(topo)
            cur = sum(times[l] for l in leaves) / len(leaves) if leaves else 0.0
            if cur > 0:
                s = H / cur
                topo = TreeState(
                    node_ids=topo.node_ids,
                    root_id=topo.root_id,
                    edges=topo.edges,
                    branch_lengths={e: v * s for e, v in topo.branch_lengths.items()},
                    node_seqs={},
                    active_leaves=leaves,
                )

        seqs = evolve_pyvolve(topo, root_seq, model=self.seq_model, seed=seed)
        seqs[topo.root_id] = root_seq
        # Table 5 metrics use leaf pool only — drop internals from payload.
        leaves = set(T.leaf_labels(topo))
        out = {topo.root_id: root_seq}
        for lid in leaves:
            if lid in seqs:
                out[lid] = seqs[lid]
        return out


def _induced_subtree(topo, keep_leaves: set[str]):
    """Prune TreeState to the subtree spanning ``keep_leaves`` (plus root path)."""
    from collections import defaultdict, deque

    from benchmarks.metrics import trees as T
    from src.tree_state import TreeState

    children = T.children_map(topo)
    parent = {c: p for p, c in topo.edges}

    # nodes on paths from kept leaves to root
    needed = set(keep_leaves)
    for leaf in keep_leaves:
        n = leaf
        while n != topo.root_id and n in parent:
            n = parent[n]
            needed.add(n)
        needed.add(topo.root_id)

    # rebuild edges among needed nodes, contracting degree-2 chains except root
    # First: keep only edges where both ends needed
    raw_children: dict[str, list[str]] = defaultdict(list)
    for p, c in topo.edges:
        if p in needed and c in needed:
            raw_children[p].append(c)

    # Contract unary internal nodes (not root, not kept leaf)
    def is_kept_leaf(n: str) -> bool:
        return n in keep_leaves

    new_edges: list[tuple[str, str]] = []
    new_bls: dict[tuple[str, str], float] = {}

    def walk(u: str, bl_acc: float, from_parent: str | None):
        kids = raw_children.get(u, [])
        if u == topo.root_id:
            for v in kids:
                bl = float(topo.branch_lengths.get((u, v), 0.0))
                walk(v, bl, u)
            return
        # if unary and not a kept leaf: skip node, accumulate BL
        if len(kids) == 1 and not is_kept_leaf(u):
            v = kids[0]
            bl = bl_acc + float(topo.branch_lengths.get((u, v), 0.0))
            walk(v, bl, from_parent)
            return
        # emit edge from_parent → u
        if from_parent is not None:
            new_edges.append((from_parent, u))
            new_bls[(from_parent, u)] = bl_acc
        for v in kids:
            bl = float(topo.branch_lengths.get((u, v), 0.0))
            walk(v, bl, u)

    walk(topo.root_id, 0.0, None)
    if not new_edges:
        return None

    # collect nodes
    nodes = {topo.root_id}
    for p, c in new_edges:
        nodes.add(p)
        nodes.add(c)
    leaves = [n for n in nodes if n in keep_leaves]
    if len(leaves) != len(keep_leaves):
        return None
    # BFS order node_ids
    node_ids, seen, q = [topo.root_id], {topo.root_id}, deque([topo.root_id])
    cm: dict[str, list[str]] = defaultdict(list)
    for p, c in new_edges:
        cm[p].append(c)
    while q:
        for c in cm.get(q.popleft(), []):
            if c not in seen:
                seen.add(c)
                node_ids.append(c)
                q.append(c)
    return TreeState(
        node_ids=node_ids,
        root_id=topo.root_id,
        edges=new_edges,
        branch_lengths=new_bls,
        node_seqs={},
        active_leaves=leaves,
    )
