"""Unit tests for benchmarks/heldout/build_examples.py (induced subtree + collapse)."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.tree_state import TreeState
from benchmarks.heldout import build_examples as B
from benchmarks.metrics import trees as T


def _tree():
    #   r -(1)-> a -(1)-> b -(1)-> d(leaf)
    #                       -(1)-> e(leaf)
    #   r -(2)-> c(leaf)
    # a is a pass-through (single child b) -> should collapse.
    ids = ["r", "a", "b", "c", "d", "e"]
    edges = [("r", "a"), ("a", "b"), ("b", "d"), ("b", "e"), ("r", "c")]
    bls = {("r", "a"): 1.0, ("a", "b"): 1.0, ("b", "d"): 1.0,
           ("b", "e"): 1.0, ("r", "c"): 2.0}
    seqs = {n: f"AA{n}A" for n in ids}
    return TreeState(node_ids=ids, root_id="r", edges=edges,
                     branch_lengths=bls, node_seqs=seqs,
                     active_leaves=["c", "d", "e"])


def test_collapse_all_leaves():
    sub = B.induced_subtree(_tree(), "r", ["d", "e", "c"])
    assert set(T.leaf_labels(sub)) == {"d", "e", "c"}
    e = set(sub.edges)
    assert e == {("r", "b"), ("r", "c"), ("b", "d"), ("b", "e")}
    # a collapsed into r->b: bl 1+1 = 2
    assert abs(sub.branch_lengths[("r", "b")] - 2.0) < 1e-9
    assert abs(sub.branch_lengths[("r", "c")] - 2.0) < 1e-9
    assert abs(sub.branch_lengths[("b", "d")] - 1.0) < 1e-9


def test_collapse_chain_two_leaves():
    sub = B.induced_subtree(_tree(), "r", ["d", "c"])
    assert set(T.leaf_labels(sub)) == {"d", "c"}
    assert set(sub.edges) == {("r", "d"), ("r", "c")}
    # a and b both collapse into r->d: bl 1+1+1 = 3
    assert abs(sub.branch_lengths[("r", "d")] - 3.0) < 1e-9
    assert abs(sub.branch_lengths[("r", "c")] - 2.0) < 1e-9


def test_root_preserved_and_seqs_carried():
    sub = B.induced_subtree(_tree(), "r", ["d", "c"])
    assert sub.root_id == "r"
    assert sub.node_seqs["r"] == "AArA"       # carried from original
    assert sub.node_seqs["d"] == "AAdA"


def test_root_to_tip_and_H():
    sub = B.induced_subtree(_tree(), "r", ["d", "c"])
    rtt = sorted(B.root_to_tip(sub))
    assert rtt == [2.0, 3.0]                   # c at 2, d at 3
    assert abs(sum(rtt) / 2 - 2.5) < 1e-9      # H = mean root-to-tip


def test_exactly_N_leaves():
    # subselecting 2 of 3 leaves yields exactly 2 terminal leaves
    sub = B.induced_subtree(_tree(), "r", ["e", "c"])
    assert len(T.leaf_labels(sub)) == 2


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    failed = 0
    for fn in fns:
        try:
            fn(); print(f"PASS  {fn.__name__}")
        except Exception as e:  # noqa: BLE001
            failed += 1; print(f"FAIL  {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
