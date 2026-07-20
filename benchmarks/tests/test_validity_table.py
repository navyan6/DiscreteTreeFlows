"""Unit tests for benchmarks/validity.py and benchmarks/make_table.py."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np

from src.tree_state import TreeState
from benchmarks import validity as V
from benchmarks import make_table as MT


def _valid_tree(root_seq="MMMM", bl=0.05):
    ids = ["R", "l1", "l2"]
    return TreeState(node_ids=ids, root_id="R",
                     edges=[("R", "l1"), ("R", "l2")],
                     branch_lengths={("R", "l1"): bl, ("R", "l2"): bl},
                     node_seqs={"R": root_seq, "l1": "MMMM", "l2": "WWWW"},
                     active_leaves=["l1", "l2"])


def test_valid_passes():
    r = V.validate(_valid_tree(), root_seq="MMMM", N=2, H=0.05)
    assert r["valid"], r["reasons"]
    assert r["n_leaves"] == 2
    assert abs(r["mean_root_to_tip"] - 0.05) < 1e-9


def test_wrong_N_flagged():
    r = V.validate(_valid_tree(), root_seq="MMMM", N=3, H=0.05)
    assert not r["valid"]
    assert any("n_leaves" in x for x in r["reasons"])


def test_root_seq_mismatch_flagged():
    r = V.validate(_valid_tree(root_seq="MMMM"), root_seq="WWWW", N=2, H=0.05)
    assert not r["valid"] and "root_seq_mismatch" in r["reasons"]


def test_horizon_tol_flagged():
    r = V.validate(_valid_tree(bl=0.05), root_seq="MMMM", N=2, H=1.0)  # asked H=1, got 0.05
    assert not r["valid"] and "horizon_out_of_tol" in r["reasons"]


def test_bad_alphabet_flagged():
    t = _valid_tree()
    t.node_seqs["l1"] = "1234"
    r = V.validate(t, root_seq="MMMM", N=2, H=0.05)
    assert not r["valid"] and "bad_alphabet" in r["reasons"]


# ── make_table ───────────────────────────────────────────────────────────────

def _results_rows():
    rows = []
    for method, rf in [("treesbm", 0.30), ("neutral_bd", 0.50)]:
        for root in range(5):
            rows.append(dict(method=method, track="empirical", N=16, root_id=root,
                             valid=5, tree_kl="nan", split_kl="nan", rf=rf,
                             quartet=0.20, branch_w_all=0.10, terminal_edit=0.05))
    return rows


def test_summarize_means_and_se():
    s = MT.summarize(_results_rows(), ci="se")
    assert {r["method"] for r in s} == {"treesbm", "neutral_bd"}
    ts = next(r for r in s if r["method"] == "treesbm")
    assert abs(ts["rf"] - 0.30) < 1e-9
    assert abs(ts["rf_unc"]) < 1e-9          # identical values -> SE 0
    assert ts["n_roots"] == 5


def test_latex_contains_methods():
    tex = MT.to_latex(MT.summarize(_results_rows(), ci="se"), "test")
    assert "treesbm" in tex and r"\toprule" in tex and "Tree-KL" in tex


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
