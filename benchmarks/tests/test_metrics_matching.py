"""Unit tests for benchmarks/metrics/matching.py."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np

from src.tree_state import TreeState
from benchmarks.metrics import matching as M


def _tree(seqs):
    # root -> a,b ; a -> c,d  (leaves b,c,d), unit branch lengths
    ids = ["root", "a", "b", "c", "d"]
    return TreeState(
        node_ids=ids, root_id="root",
        edges=[("root", "a"), ("root", "b"), ("a", "c"), ("a", "d")],
        branch_lengths={("root", "a"): 1.0, ("root", "b"): 1.0,
                        ("a", "c"): 1.0, ("a", "d"): 1.0},
        node_seqs=seqs, active_leaves=["b", "c", "d"],
    )


def test_patristic_matrix():
    t = _tree({n: "AA" for n in ["root", "a", "b", "c", "d"]})
    leaves, Pm = M.patristic_matrix(t)
    idx = {l: i for i, l in enumerate(leaves)}
    # c,d share parent a: patristic = 1+1 = 2 ; b,c: via root = time(b)=1 + time(c)=2 - 2*time(root)=0 => 3
    assert abs(Pm[idx["c"], idx["d"]] - 2.0) < 1e-9
    assert abs(Pm[idx["b"], idx["c"]] - 3.0) < 1e-9


def test_hungarian_match_recovers_permutation():
    r = M.hungarian_match(true_seqs=["AAA", "BBB"], gen_seqs=["BBB", "AAA"])
    assert r["mean_matched_identity"] == 1.0
    assert r["mean_matched_hamming"] == 0.0
    assert r["n_matched"] == 2


def test_seq_patristic_correlation_positive_when_coupled():
    # sequences whose divergence tracks tree distance: c,d close; b far
    seqs = {"root": "AAAAAA", "a": "AAAAAA",
            "b": "TTTTTT",       # far from c,d
            "c": "AAAAAA", "d": "AAAAAT"}  # c,d nearly identical (close)
    t = _tree(seqs)
    corr = M.seq_patristic_correlation(t)
    assert corr > 0.5  # closer in tree => closer in sequence


def test_matched_patristic_agreement_identical_trees():
    seqs = {"root": "AAAA", "a": "AAAA", "b": "TTTT", "c": "AAAA", "d": "AACA"}
    t = _tree(seqs)
    res = M.matched_patristic_agreement(t, t)
    # matching a tree to itself: perfect patristic correlation, ~0 nrmse
    assert res["patristic_corr"] > 0.99
    assert res["patristic_nrmse"] < 1e-6


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
