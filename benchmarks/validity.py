"""
Validity checks for every generated sample (Validity checks section).

A sample is only scored if it is a well-formed rooted tree matching the task
request. Failures are recorded (never silently dropped) so the table can report
attempted vs valid counts and failure reasons.
"""

from __future__ import annotations

import math
from collections import Counter, deque

from src.tree_state import TreeState
from benchmarks.metrics import trees as T

AA_VOCAB = set("ACDEFGHIKLMNPQRSTVWY")
ALLOWED = AA_VOCAB | {"-"}


def validate(gen: TreeState, root_seq: str, N: int, H: float,
             h_rel_tol: float = 0.5) -> dict:
    """
    Returns {valid, reasons, n_leaves, mean_root_to_tip}. `h_rel_tol` is the
    allowed relative deviation of mean root-to-tip distance from the requested H.
    """
    reasons: list[str] = []
    L = len(root_seq)

    # duplicate node ids
    if len(gen.node_ids) != len(set(gen.node_ids)):
        reasons.append("duplicate_node_ids")

    # exactly one root (node with no parent), and it is gen.root_id
    parents = {c for _, c in gen.edges}
    roots = [n for n in gen.node_ids if n not in parents]
    if roots != [gen.root_id]:
        reasons.append(f"root_set={roots[:3]}")

    # every non-root node has exactly one parent
    pc = Counter(c for _, c in gen.edges)
    if any(pc.get(n, 0) != 1 for n in gen.node_ids if n != gen.root_id):
        reasons.append("non_root_parent_count_not_1")

    # connected acyclic tree: |E| == |V|-1 and all nodes reachable from root
    if len(gen.edges) != len(gen.node_ids) - 1:
        reasons.append("edge_count_not_tree")
    cm = T.children_map(gen)
    seen, q = {gen.root_id}, deque([gen.root_id])
    while q:
        for c in cm.get(q.popleft(), []):
            if c not in seen:
                seen.add(c); q.append(c)
    if len(seen) != len(gen.node_ids):
        reasons.append("disconnected")

    # exactly N terminal leaves
    leaves = T.leaf_labels(gen)
    if len(leaves) != N:
        reasons.append(f"n_leaves_{len(leaves)}_ne_{N}")

    # branch lengths nonnegative + finite
    if any((not math.isfinite(v)) or v < 0 for v in gen.branch_lengths.values()):
        reasons.append("bad_branch_length")

    # root sequence matches the supplied root exactly
    if gen.node_seqs.get(gen.root_id) != root_seq:
        reasons.append("root_seq_mismatch")

    # sequence length + alphabet
    if any(len(s) != L for s in gen.node_seqs.values()):
        reasons.append("seq_length")
    if any(set(s) - ALLOWED for s in gen.node_seqs.values()):
        reasons.append("bad_alphabet")

    # horizon within relative tolerance
    times = T.node_times(gen)
    rtt = [times[l] for l in leaves]
    mean_rtt = sum(rtt) / len(rtt) if rtt else 0.0
    if H > 0 and abs(mean_rtt - H) > h_rel_tol * H:
        reasons.append("horizon_out_of_tol")

    return {"valid": len(reasons) == 0, "reasons": reasons,
            "n_leaves": len(leaves), "mean_root_to_tip": mean_rtt}
