"""
Ground-truth extraction for Track A tasks, straight from a real tree — no
simulator. A `TreeState` plus per-node times (branch-length cumulative, or real
`numdate`s) fully determines the event history we forecast against.

Provides:
  event_stream      chronological list of branch / mutation / stop events
  lineages_at       active lineages (edges crossing time t)
  next_event_after  the true next event after time t  (A3 label)
  reveal_prefix     deterministic partial tree revealed up to a time fraction,
                    with in-flight lineages cut into tips  (A2 input + remainder)

Mutation timing granularity: parent/child sequences give *which* sites mutate on
an edge, not *when*, so edge mutations are attributed to the child-node event.
"""

from __future__ import annotations

from src.tree_state import TreeState
from benchmarks.metrics import trees as T

__all__ = [
    "event_stream", "lineages_at", "next_event_after", "reveal_prefix",
]


def _edge_mutations(seq_p: str, seq_c: str) -> list[tuple[int, str, str]]:
    L = min(len(seq_p), len(seq_c))
    return [(i, seq_p[i], seq_c[i]) for i in range(L) if seq_p[i] != seq_c[i]]


def event_stream(tree: TreeState, node_times: dict[str, float] | None = None) -> list[dict]:
    """
    Chronological events. For each edge (p, c) at time(c): a `branch` (c internal)
    or `stop` (c leaf) event on lineage p, carrying the edge's mutations.
    """
    times = node_times if node_times is not None else T.node_times(tree)
    cm = T.children_map(tree)
    events: list[dict] = []
    for p, c in tree.edges:
        muts = _edge_mutations(tree.node_seqs[p], tree.node_seqs[c])
        events.append({
            "time": times[c],
            "waiting_from_parent": times[c] - times[p],
            "lineage": p,
            "child": c,
            "type": "branch" if cm.get(c) else "stop",
            "mutations": muts,                       # list of (pos, from, to)
            "sites": {m[0] for m in muts},
        })
    events.sort(key=lambda e: (e["time"], e["child"]))
    return events


def lineages_at(tree: TreeState, t: float, node_times: dict[str, float] | None = None
                ) -> list[tuple[str, str]]:
    """Edges (p, c) crossing time t: time(p) <= t < time(c) — the active lineages."""
    times = node_times if node_times is not None else T.node_times(tree)
    return [(p, c) for p, c in tree.edges if times[p] <= t < times[c]]


def next_event_after(tree: TreeState, t: float, node_times: dict[str, float] | None = None
                     ) -> dict | None:
    """
    The true next event strictly after time t: among lineages crossing t, the edge
    whose child has the smallest time. Returns None if t is past all events.
    (A3 ground-truth label at an intermediate state.)
    """
    times = node_times if node_times is not None else T.node_times(tree)
    crossing = lineages_at(tree, t, times)
    if not crossing:
        return None
    p, c = min(crossing, key=lambda e: (times[e[1]], e[1]))
    muts = _edge_mutations(tree.node_seqs[p], tree.node_seqs[c])
    cm = T.children_map(tree)
    return {
        "time": times[c],
        "waiting_time": times[c] - t,
        "lineage": p,                                # which active lineage
        "child": c,
        "type": "branch" if cm.get(c) else "stop",
        "mutations": muts,
        "sites": {m[0] for m in muts},
        "n_active": len(crossing),
    }


def reveal_prefix(tree: TreeState, fraction: float,
                  node_times: dict[str, float] | None = None
                  ) -> dict:
    """
    Reveal the chronological prefix up to cutoff = t_root + fraction*(t_max-t_root).

    Nodes with time <= cutoff are kept with their TRUE sequences; each edge (p, c)
    crossing the cutoff is cut into a tip `f"{c}__cut"` on lineage p (carrying p's
    sequence) — these tips are the active lineages the model must continue.

    Returns:
      partial      TreeState (revealed nodes + cut tips; active_leaves = cut tips)
      cutoff       the time cutoff
      future_nodes list of hidden node ids (time > cutoff)
      remainder_by_lineage  {cut_tip_id: [true future node ids on that lineage]}
    """
    times = node_times if node_times is not None else T.node_times(tree)
    t_root = times[tree.root_id]
    t_max = max(times.values())
    cutoff = t_root if t_max <= t_root else t_root + fraction * (t_max - t_root)

    retained = {n for n in tree.node_ids if times[n] <= cutoff}
    retained.add(tree.root_id)

    node_ids = [n for n in tree.node_ids if n in retained]
    edges = [(p, c) for p, c in tree.edges if p in retained and c in retained]
    branch_lengths = {(p, c): tree.branch_lengths.get((p, c), 0.0) for p, c in edges}
    seqs = {n: tree.node_seqs[n] for n in node_ids}

    # cut in-flight lineages into tips
    cm_full = T.children_map(tree)
    crossing = [(p, c) for p, c in tree.edges if p in retained and c not in retained]
    active: list[str] = []
    remainder: dict[str, list[str]] = {}
    for p, c in crossing:
        tip = f"{c}__cut"
        node_ids.append(tip)
        edges.append((p, tip))
        branch_lengths[(p, tip)] = max(0.0, cutoff - times[p])
        seqs[tip] = tree.node_seqs[p]            # sequence at cutoff (edge muts attributed to child)
        active.append(tip)
        # true future subtree hanging off this cut lineage
        sub, stack = [], [c]
        while stack:
            n = stack.pop()
            sub.append(n)
            stack.extend(cm_full.get(n, []))
        remainder[tip] = sub

    partial = TreeState(
        node_ids=node_ids, root_id=tree.root_id, edges=edges,
        branch_lengths=branch_lengths, node_seqs=seqs,
        active_leaves=active if active else [tree.root_id],
    )
    future_nodes = [n for n in tree.node_ids if n not in retained]
    return {
        "partial": partial,
        "cutoff": cutoff,
        "future_nodes": future_nodes,
        "remainder_by_lineage": remainder,
    }
