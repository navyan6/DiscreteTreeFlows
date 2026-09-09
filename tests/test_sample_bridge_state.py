"""sample_bridge_state must return a root-connected tree."""

from src.bridge.sample_bridge_state import sample_bridge_state
from src.tree_state import TreeState
from src.treeencoder.structural_features import compute_depth


def test_bridge_drops_unconnected_zero_time_descendants():
    """Child with fallback numdate=0 below a later parent must not stay in T_t."""
    node_ids = ["root", "a", "b"]
    times = {"root": 0.0, "a": 1.0, "b": 0.0}
    edges = [("root", "a"), ("a", "b")]
    bls = {("root", "a"): 1.0, ("a", "b"): -1.0}
    seqs = {n: "AA" for n in node_ids}
    T = sample_bridge_state(
        t=0.0,
        node_ids=node_ids,
        node_times_dict=times,
        edges=edges,
        branch_lengths=bls,
        seqs=seqs,
        root_id="root",
    )
    assert T["node_ids_t"] == ["root"]
    tree = TreeState(
        node_ids=T["node_ids_t"],
        root_id="root",
        edges=T["edges_t"],
        branch_lengths=T["branch_lengths_t"],
        node_seqs=T["seqs_t"],
        active_leaves=T["active_leaves_t"],
    )
    assert compute_depth(tree) == {"root": 0}


def test_bridge_full_tree_at_t1():
    node_ids = ["root", "a", "b"]
    times = {"root": 0.0, "a": 0.5, "b": 1.0}
    edges = [("root", "a"), ("a", "b")]
    bls = {("root", "a"): 0.5, ("a", "b"): 0.5}
    seqs = {n: "AA" for n in node_ids}
    T = sample_bridge_state(
        t=1.0,
        node_ids=node_ids,
        node_times_dict=times,
        edges=edges,
        branch_lengths=bls,
        seqs=seqs,
        root_id="root",
    )
    assert set(T["node_ids_t"]) == {"root", "a", "b"}
    tree = TreeState(
        node_ids=T["node_ids_t"],
        root_id="root",
        edges=T["edges_t"],
        branch_lengths=T["branch_lengths_t"],
        node_seqs=T["seqs_t"],
        active_leaves=T["active_leaves_t"],
    )
    depths = compute_depth(tree)
    assert depths == {"root": 0, "a": 1, "b": 2}
