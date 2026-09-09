"""Unit tests for matched.py, branch_lengths.py, and JS variants in distributions.py."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.tree_state import TreeState
from benchmarks.metrics import matched as M
from benchmarks.metrics import branch_lengths as BL
from benchmarks.metrics import distributions as D


def _quartet(a, b, c, d, seqs, bl=1.0, ids=("R", "A", "B")):
    # rooted 4-leaf tree: (a,b),(c,d)
    R, IA, IB = ids
    node_ids = [R, IA, IB, a, b, c, d]
    edges = [(R, IA), (R, IB), (IA, a), (IA, b), (IB, c), (IB, d)]
    bls = {e: bl for e in edges}
    ns = {R: "----", IA: "----", IB: "----", a: seqs[0], b: seqs[1], c: seqs[2], d: seqs[3]}
    return TreeState(node_ids=node_ids, root_id=R, edges=edges,
                     branch_lengths=bls, node_seqs=ns, active_leaves=[a, b, c, d])


REF = _quartet("l1", "l2", "l3", "l4", ["MMMM", "WWWW", "YYYY", "FFFF"])


def test_matched_rf_zero_when_isomorphic():
    # gen leaves g1..g4 carry the same seqs as l1..l4, same topology -> RF 0
    gen = _quartet("g1", "g2", "g3", "g4", ["MMMM", "WWWW", "YYYY", "FFFF"],
                   ids=("r", "a", "b"))
    assert M.sequence_matched_rf(gen, REF) == 0.0


def test_matched_rf_positive_on_different_topology():
    # gen groups {l1,l3} and {l2,l4} instead of {l1,l2},{l3,l4}
    gen = _quartet("g1", "g3", "g2", "g4", ["MMMM", "YYYY", "WWWW", "FFFF"],
                   ids=("r", "a", "b"))
    rf = M.sequence_matched_rf(gen, REF)
    assert rf > 0.0


def test_terminal_edit_distance():
    gen = _quartet("g1", "g2", "g3", "g4", ["MMMM", "WWWW", "YYYY", "FFFA"],
                   ids=("r", "a", "b"))   # g4 differs from l4 by 1/4
    res = M.terminal_edit_distance(gen, REF)
    assert abs(res["mean"] - (0.25 / 4)) < 1e-9   # only one leaf off by 1 of 4


def test_branch_length_wasserstein():
    gen = _quartet("g1", "g2", "g3", "g4", ["MMMM", "WWWW", "YYYY", "FFFF"],
                   bl=0.1, ids=("r", "a", "b"))
    w = BL.branch_length_wasserstein(gen, REF)   # ref bl=1.0, gen bl=0.1
    assert abs(w["all"] - 0.9) < 1e-9


def test_tree_and_split_js():
    same = [REF for _ in range(4)]
    assert D.tree_js(same, same) < 1e-9
    assert D.split_js(same, same) < 1e-6
    other = _quartet("l1", "l3", "l2", "l4", ["MMMM", "YYYY", "WWWW", "FFFF"])
    assert D.tree_js([REF] * 4, [other] * 4) > 0.0
    assert D.split_js([REF] * 4, [other] * 4) > 0.0


def test_quartet_guarded():
    # tqdist may not be installed locally; either a float or a clean ImportError
    gen = _quartet("g1", "g2", "g3", "g4", ["MMMM", "WWWW", "YYYY", "FFFF"],
                   ids=("r", "a", "b"))
    try:
        q = M.quartet_distance(gen, REF)
        assert isinstance(q, float)
    except ImportError:
        pass


def _stem(tree: TreeState, stem_id="STEM", bl=0.01) -> TreeState:
    """Wrap `tree` with a new unary root -> stem_id -> old root (dendropy-style)."""
    return TreeState(
        node_ids=[stem_id] + tree.node_ids, root_id=stem_id,
        edges=[(stem_id, tree.root_id)] + tree.edges,
        branch_lengths={(stem_id, tree.root_id): bl, **tree.branch_lengths},
        node_seqs=tree.node_seqs, active_leaves=list(tree.active_leaves),
    )


def test_drop_unary_root_noop_when_not_unary():
    assert M._drop_unary_root(REF).root_id == REF.root_id


def test_drop_unary_root_promotes_child():
    stemmed = _stem(REF)
    dropped = M._drop_unary_root(stemmed)
    assert dropped.root_id == REF.root_id
    assert set(dropped.node_ids) == set(REF.node_ids)
    assert "STEM" not in dropped.node_ids


def test_quartet_unary_root_matches_isomorphic():
    # gen has a dendropy-style stem root; topology is otherwise isomorphic to
    # REF, so quartet distance should be 0 (not a "leaves don't agree" abort).
    gen = _stem(_quartet("g1", "g2", "g3", "g4", ["MMMM", "WWWW", "YYYY", "FFFF"],
                        ids=("r", "a", "b")))
    try:
        q = M.quartet_distance(gen, REF)
        assert q == 0.0
    except ImportError:
        pass


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
