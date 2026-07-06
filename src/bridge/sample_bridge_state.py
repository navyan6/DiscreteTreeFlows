"""
Algorithm 2: SampleBridgeState

Interpolates between T0 = {x0} (root-only) and T1 (observed tree) at time t in [0,1].

Each branch (u,v) is retained if time(u) <= t_cut = t_root + t*(t_max - t_root).
Branch lengths are scaled by t.
Mutations along each branch are Bernoulli-sampled with probability proportional
to how far t_cut has advanced along the branch's time interval.
"""

import random


def sample_bridge_state(
    t: float,
    node_ids: list[str],
    node_times_dict: dict[str, float],
    edges: list[tuple[str, str]],
    branch_lengths: dict[tuple[str, str], float],
    seqs: dict[str, str],
    root_id: str,
) -> dict:
#    Sample a bridge state T_t at time t in [0,1] between T0 and T1.
    # 1. Time cutoff
    times = [node_times_dict[nid] for nid in node_ids]
    t_root_val = min(times)
    t_max_val = max(times)

    t_cut = t_root_val if t_max_val <= t_root_val else (
        t_root_val + t * (t_max_val - t_root_val)
    )

    # 2. Retain nodes whose numdate <= t_cut (always keep root)
    retained = {nid for nid in node_ids if node_times_dict[nid] <= t_cut}
    retained.add(root_id)

    # 3. Retain edges where both parent and child are retained
    retained_edges = [(p, c) for p, c in edges if p in retained and c in retained]

    # 4. Branch lengths for retained edges are their full T1 values (branch is fully elapsed)
    branch_lengths_t = {(p, c): branch_lengths[(p, c)] for p, c in retained_edges}

    # 5. Partial sequences via Bernoulli-sampling mutations per branch
    parent_map = {c: p for p, c in retained_edges}
    node_ids_t = [nid for nid in node_ids if nid in retained]  # preserve BFS order

    seqs_t = {root_id: seqs[root_id]}
    for nid in node_ids_t:
        if nid == root_id:
            continue
        parent = parent_map.get(nid)
        if parent is None or parent not in seqs_t:
            seqs_t[nid] = seqs[root_id]
            continue

        p_seq = seqs_t[parent]
        c_seq = seqs[nid]

        # Fraction of branch time elapsed at t_cut
        t_p = node_times_dict[parent]
        t_c = node_times_dict[nid]
        frac = min(1.0, max(0.0, (t_cut - t_p) / (t_c - t_p))) if t_c > t_p else 1.0

        # Each differing position mutated independently with prob frac
        partial = list(p_seq)
        L = min(len(p_seq), len(c_seq))
        for pos in range(L):
            if p_seq[pos] != c_seq[pos] and random.random() < frac:
                partial[pos] = c_seq[pos]
        seqs_t[nid] = "".join(partial)

    # 6. Active leaves of T_t 
    has_children_t = {p for p, c in retained_edges}
    active_leaves_t = [nid for nid in node_ids_t if nid not in has_children_t]

    # 7. T1 child info for each active leaf (supervision targets)
    T1_children_map: dict[str, list[str]] = {}
    for p, c in edges:
        T1_children_map.setdefault(p, []).append(c)

    T1_child_counts = {v: len(T1_children_map.get(v, [])) for v in active_leaves_t}
    T1_child_bls = {
        v: [branch_lengths[(v, c)] for c in T1_children_map.get(v, [])]
        for v in active_leaves_t
    }

    return {
        "node_ids_t": node_ids_t,
        "edges_t": retained_edges,
        "branch_lengths_t": branch_lengths_t,
        "seqs_t": seqs_t,
        "active_leaves_t": active_leaves_t,
        "T1_child_counts": T1_child_counts,
        "T1_child_bls": T1_child_bls,
        "t_cut": t_cut,
    }
