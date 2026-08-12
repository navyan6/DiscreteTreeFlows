"""Pure-Python checks for coverage-curve helpers (no dendropy/ESM)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.tree_state import TreeState
from benchmarks.coverage_curves import (
    k_grid, score_pool, emerging_clade_fingerprints, clade_recall_at_k,
)
from benchmarks.metrics import sequences as S


def test_k_grid():
    assert k_grid(100, 10) == list(range(10, 101, 10))
    assert k_grid(25, 10) == [10, 20, 25]


def test_score_pool_monotone_edit_and_coverage():
    root = "AAAA"
    gt = ["ABAA", "AABA"]
    # exact cover of both targets
    pool_full = ["ABAA", "AABA", "AAAA"]
    m = score_pool(gt, root, pool_full, eps_frac=0.0, max_gt=10, max_gen=10, seed=0,
                   e_list=[0, 1, 2])
    assert m["coverage"] == 1.0
    assert m["coverage_obs_e0"] == 1.0
    assert m["mean_min_edit"] == 0.0
    assert m["site_recall"] == 1.0  # sites {1,2} recovered

    # worse pool: only one close sequence
    pool_partial = ["ABAA", "AAAA"]
    m2 = score_pool(gt, root, pool_partial, eps_frac=0.0, max_gt=10, max_gen=10, seed=0,
                    e_list=[0, 1, 2])
    assert m2["coverage"] == 0.5
    assert m2["coverage_obs_e0"] == 0.5
    assert m2["coverage_obs_e1"] == 1.0  # second tip is 1 mut from AAAA? AABA vs ABAA = 2; vs AAAA = 1
    assert m2["mean_min_edit"] > m["mean_min_edit"]
    # gen→obs: both gen seqs are within e=0 of some obs (ABAA exact; AAAA within 1 of both)
    assert m2["frac_gen_e0"] == 0.5  # only ABAA exact-matches an obs tip


def test_coverage_at_e_and_frac_gen():
    assert S.coverage_at_e(["AAA"], ["AAA"], e=0) == 1.0
    assert S.coverage_at_e(["AAA"], ["ABA"], e=0) == 0.0
    assert S.coverage_at_e(["AAA"], ["ABA"], e=1) == 1.0
    assert S.frac_gen_within_e(["AAA"], ["AAA", "ABA"], e=0) == 0.5
    assert S.frac_gen_within_e(["AAA"], ["AAA", "ABA"], e=1) == 1.0


def test_clade_recall_fingerprint():
    # Root AAAA; cherry sharing mut at pos 0 -> B***; other leaf private mut
    # Tree: root -> (i1 -> (L1, L2), L3)
    # L1=BCAA, L2=BADA share B at 0; L3=AAAE
    root = "AAAA"
    edges = [("root", "i1"), ("i1", "L1"), ("i1", "L2"), ("root", "L3")]
    bl = {e: 1.0 for e in edges}
    seqs = {"root": root, "i1": "BAAA", "L1": "BCAA", "L2": "BADA", "L3": "AAAE"}
    gt = TreeState(
        node_ids=list(seqs),
        root_id="root",
        edges=edges,
        branch_lengths=bl,
        node_seqs=seqs,
        active_leaves=["L1", "L2", "L3"],
    )
    fps = emerging_clade_fingerprints(gt, root, min_clade_size=2, min_shared_muts=1)
    assert any((0, "B") in fp for fp in fps)

    # Gen with shared B mut covers the clade
    cr, n = clade_recall_at_k(fps, ["BZZZ"], root)
    assert n >= 1
    assert cr == 1.0

    # Gen without shared fingerprint misses
    cr0, _ = clade_recall_at_k(fps, ["AAAZ"], root)
    assert cr0 == 0.0


if __name__ == "__main__":
    test_k_grid()
    test_score_pool_monotone_edit_and_coverage()
    test_coverage_at_e_and_frac_gen()
    test_clade_recall_fingerprint()
    print("ok")
