"""Unit tests for benchmarks/metrics/events.py and benchmarks/task_data.py."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np

from src.tree_state import TreeState
from benchmarks.metrics import events as E
from benchmarks import task_data as TD


# ── events.py ───────────────────────────────────────────────────────────────

def test_event_nll_and_brier_perfect():
    probs = np.array([[1.0, 0.0], [0.0, 1.0]])
    true = np.array([0, 1])
    assert E.event_nll(probs, true) < 1e-6
    assert E.brier_score(probs, true) < 1e-9


def test_event_nll_uniform():
    probs = np.array([[0.5, 0.5], [0.5, 0.5]])
    true = np.array([0, 1])
    assert abs(E.event_nll(probs, true) - np.log(2)) < 1e-9


def test_ece_perfectly_calibrated_is_zero():
    # confidence 1.0 and always correct -> ECE 0
    conf = np.array([1.0, 1.0, 1.0])
    correct = np.array([1, 1, 1])
    assert E.expected_calibration_error(conf, correct) < 1e-9
    # overconfident: conf 1.0 but half wrong -> ECE 0.5
    conf2 = np.array([1.0, 1.0, 1.0, 1.0])
    correct2 = np.array([1, 0, 1, 0])
    assert abs(E.expected_calibration_error(conf2, correct2) - 0.5) < 1e-9


def test_time_to_event_error():
    r = E.time_to_event_error([1.0, 2.0], [1.5, 2.0])
    assert abs(r["mae"] - 0.25) < 1e-9


def test_mutation_site_recall_at_k():
    scores = np.array([0.1, 0.9, 0.2, 0.8, 0.05])  # top sites: 1, 3, ...
    r = E.mutation_site_recall_at_k(scores, true_sites={1, 4}, ks=(1, 2, 5))
    assert r[1] == 0.5     # top1={1}; recovers 1 of {1,4}
    assert r[2] == 0.5     # top2={1,3}; still just 1
    assert r[5] == 1.0     # all sites -> both recovered


# ── task_data.py ────────────────────────────────────────────────────────────

def _timed_tree():
    # times via branch lengths: root0 a1 c2 d4 b5
    #   root -(1)-> a -(1)-> c(leaf)
    #             a -(3)-> d(leaf)   [mut pos2 A->C]
    #   root -(5)-> b(leaf)
    #   root->a edge mutates pos1 A->B
    ids = ["root", "a", "b", "c", "d"]
    seqs = {"root": "AAAA", "a": "ABAA", "b": "AAAA", "c": "ABAA", "d": "ABCA"}
    return TreeState(
        node_ids=ids, root_id="root",
        edges=[("root", "a"), ("root", "b"), ("a", "c"), ("a", "d")],
        branch_lengths={("root", "a"): 1.0, ("root", "b"): 5.0,
                        ("a", "c"): 1.0, ("a", "d"): 3.0},
        node_seqs=seqs, active_leaves=["b", "c", "d"],
    )


def test_event_stream_sorted_with_mutations():
    t = _timed_tree()
    ev = TD.event_stream(t)
    times = [e["time"] for e in ev]
    assert times == sorted(times)
    # root->a carries the pos1 mutation, attributed to child a at time 1
    a_ev = next(e for e in ev if e["child"] == "a")
    assert a_ev["type"] == "branch" and a_ev["sites"] == {1} and abs(a_ev["time"] - 1.0) < 1e-9
    d_ev = next(e for e in ev if e["child"] == "d")
    assert d_ev["type"] == "stop" and d_ev["sites"] == {2}


def test_lineages_and_next_event():
    t = _timed_tree()
    # at t=2.5: crossing edges root->b (child t5) and a->d (child t4)
    cross = TD.lineages_at(t, 2.5)
    assert set(cross) == {("root", "b"), ("a", "d")}
    nxt = TD.next_event_after(t, 2.5)
    assert nxt["child"] == "d"                    # d(t4) before b(t5)
    assert nxt["lineage"] == "a"
    assert abs(nxt["waiting_time"] - 1.5) < 1e-9
    assert nxt["sites"] == {2}
    assert nxt["n_active"] == 2


def test_next_event_none_past_end():
    t = _timed_tree()
    assert TD.next_event_after(t, 100.0) is None


def test_reveal_prefix_cuts_inflight_lineages():
    t = _timed_tree()
    out = TD.reveal_prefix(t, fraction=0.5)   # cutoff = 0 + 0.5*5 = 2.5
    assert abs(out["cutoff"] - 2.5) < 1e-9
    partial = out["partial"]
    # revealed nodes: root, a, c (<=2.5); cut tips for root->b and a->d
    assert set(out["remainder_by_lineage"].keys()) == {"b__cut", "d__cut"}
    assert set(partial.active_leaves) == {"b__cut", "d__cut"}
    # cut tip carries the parent sequence at cutoff
    assert partial.node_seqs["b__cut"] == t.node_seqs["root"]
    assert partial.node_seqs["d__cut"] == t.node_seqs["a"]
    # future subtree off a->d cut is just {d}; off root->b is {b}
    assert out["remainder_by_lineage"]["d__cut"] == ["d"]
    assert out["remainder_by_lineage"]["b__cut"] == ["b"]
    assert set(out["future_nodes"]) == {"b", "d"}


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"FAIL  {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
