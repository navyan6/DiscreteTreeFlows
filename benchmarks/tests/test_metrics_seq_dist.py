"""Unit tests for benchmarks/metrics/sequences.py and distributions.py."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np

from src.tree_state import TreeState
from benchmarks.metrics import sequences as S
from benchmarks.metrics import distributions as D


# ── sequences ───────────────────────────────────────────────────────────────

def test_hamming_identity():
    assert S.hamming("AAA", "ABA") == 1
    assert abs(S.identity("AAA", "ABA") - 2 / 3) < 1e-9


def test_mutations_vs_root():
    assert S.mutations_vs_root("AAA", "ABA") == {(1, "B")}
    assert S.sites_vs_root("AAA", "ABA") == {1}


def test_best_of_k_and_coverage():
    assert abs(S.best_of_k_identity("AAA", ["ABA", "AAB"]) - 2 / 3) < 1e-9
    assert S.coverage_at_k(["AAA"], ["AAA"], eps_frac=0.0) == 1.0
    assert S.coverage_at_k(["AAA"], ["ABA"], eps_frac=0.0) == 0.0
    assert S.coverage_at_k(["AAA"], ["ABA"], eps_frac=0.5) == 1.0  # 1 mismatch <= 1.5


def test_mutation_pr_f1():
    root = "AAAA"
    r = S.mutation_pr_f1(gen_seqs=["ABAA"], true_seqs=["ABAA", "AACA"], root=root)
    assert r["precision"] == 1.0
    assert abs(r["recall"] - 0.5) < 1e-9
    assert abs(r["f1"] - (2 * 1.0 * 0.5 / 1.5)) < 1e-9
    assert r["n_recovered"] == 1
    assert S.unique_mutations_recovered(["ABAA"], ["ABAA", "AACA"], root) == 1


def test_sitewise_entropy():
    ent = S.sitewise_entropy(["AA", "AB"])
    assert abs(ent[0] - 0.0) < 1e-9
    assert abs(ent[1] - np.log(2)) < 1e-9


def test_positional_recovery():
    r = S.positional_recovery("AAA", "ABA", "ABA")
    assert r["mut_recovery"] == 1.0 and r["cons_retention"] == 1.0
    assert r["mut_total"] == 1 and r["cons_total"] == 2


# ── distributions ───────────────────────────────────────────────────────────

def _balanced(bl=1.0, tag=""):
    ids = [f"root{tag}", f"a{tag}", f"b{tag}", f"c{tag}", f"d{tag}"]
    return TreeState(
        node_ids=ids, root_id=ids[0],
        edges=[(ids[0], ids[1]), (ids[0], ids[2]), (ids[1], ids[3]), (ids[1], ids[4])],
        branch_lengths={(ids[0], ids[1]): bl, (ids[0], ids[2]): bl,
                        (ids[1], ids[3]): bl, (ids[1], ids[4]): bl},
        node_seqs={n: "AA" for n in ids}, active_leaves=[ids[2], ids[3], ids[4]],
    )


def _caterpillar(bl=1.0):
    # root -> (e, x); x -> (c, d)  ... a different topology than balanced above
    return TreeState(
        node_ids=["root", "e", "x", "c", "d"], root_id="root",
        edges=[("root", "e"), ("root", "x"), ("x", "c"), ("x", "d")],
        branch_lengths={("root", "e"): bl, ("root", "x"): bl,
                        ("x", "c"): bl, ("x", "d"): bl},
        node_seqs={n: "AA" for n in ["root", "e", "x", "c", "d"]},
        active_leaves=["e", "c", "d"],
    )


def test_summarize_shape():
    v = D.summarize_tree(_balanced())
    assert v.shape == (len(D.SUMMARY_FEATURES),)


def test_identical_sets_zero_distance():
    trees = [_balanced() for _ in range(6)]
    other = [_balanced() for _ in range(6)]
    assert D.energy_distance(D.summary_matrix(trees), D.summary_matrix(other)) < 1e-9
    assert D.split_kl(trees, other) < 1e-6
    assert D.tree_kl(trees, other) < 1e-6
    w = D.wasserstein_per_feature(D.summary_matrix(trees), D.summary_matrix(other))
    assert max(w.values()) < 1e-9


def test_different_topologies_positive_split_kl():
    # same leaf labels {c,d,...} differ: balanced groups {c,d} under a; use trees
    # whose clade sets differ -> split_kl > 0
    A = [_balanced() for _ in range(5)]
    # trees where the informative clade is {b,c} instead of {c,d}
    B_ids = ["root", "a", "b", "c", "d"]
    B_tree = TreeState(
        node_ids=B_ids, root_id="root",
        edges=[("root", "a"), ("root", "d"), ("a", "b"), ("a", "c")],
        branch_lengths={("root", "a"): 1.0, ("root", "d"): 1.0,
                        ("a", "b"): 1.0, ("a", "c"): 1.0},
        node_seqs={n: "AA" for n in B_ids}, active_leaves=["b", "c", "d"],
    )
    B = [B_tree for _ in range(5)]
    assert D.split_kl(A, B) > 0.0
    assert D.tree_kl(A, B) > 0.0


def test_energy_distance_detects_shift():
    small = [_balanced(bl=0.1) for _ in range(8)]
    big = [_balanced(bl=5.0) for _ in range(8)]
    assert D.energy_distance(D.summary_matrix(small), D.summary_matrix(big)) > 0.0


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
