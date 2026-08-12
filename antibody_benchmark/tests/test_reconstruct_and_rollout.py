"""Unit tests for tree validation, reconstruction, and recursive rollout."""

from __future__ import annotations

import pandas as pd
import pytest

from antibody_benchmark.models.base import IdentityNullModel, PoissonSiteMutationNullModel
from antibody_benchmark.rollout.reconstruct_trees import reconstruct_trees_from_pcp_df
from antibody_benchmark.rollout.rollout_tree import rollout_tree
from antibody_benchmark.trees import AntibodyTree, validate_tree


def _toy_pcp_df() -> pd.DataFrame:
    # root(naive) -> i -> l1,l2 ; root -> l3 ; root -> l4
    rows = [
        dict(
            sample_id="d1",
            family="f1",
            parent_name="root",
            child_name="i",
            parent_heavy="ATGCATGCATGC",
            child_heavy="ATGTATGCATGC",
            branch_length=0.05,
            depth=1,
            distance=0.05,
            parent_is_naive=True,
            child_is_leaf=False,
            v_gene_heavy="IGHV1-1",
        ),
        dict(
            sample_id="d1",
            family="f1",
            parent_name="i",
            child_name="l1",
            parent_heavy="ATGTATGCATGC",
            child_heavy="ATGTATGAATGC",
            branch_length=0.04,
            depth=2,
            distance=0.04,
            parent_is_naive=False,
            child_is_leaf=True,
            v_gene_heavy="IGHV1-1",
        ),
        dict(
            sample_id="d1",
            family="f1",
            parent_name="i",
            child_name="l2",
            parent_heavy="ATGTATGCATGC",
            child_heavy="ATGTATGCTTGC",
            branch_length=0.03,
            depth=2,
            distance=0.03,
            parent_is_naive=False,
            child_is_leaf=True,
            v_gene_heavy="IGHV1-1",
        ),
        dict(
            sample_id="d1",
            family="f1",
            parent_name="root",
            child_name="l3",
            parent_heavy="ATGCATGCATGC",
            child_heavy="ATGCATGCATGT",
            branch_length=0.02,
            depth=1,
            distance=0.02,
            parent_is_naive=True,
            child_is_leaf=True,
            v_gene_heavy="IGHV1-1",
        ),
        dict(
            sample_id="d1",
            family="f1",
            parent_name="root",
            child_name="l4",
            parent_heavy="ATGCATGCATGC",
            child_heavy="ATACATGCATGC",
            branch_length=0.02,
            depth=1,
            distance=0.02,
            parent_is_naive=True,
            child_is_leaf=True,
            v_gene_heavy="IGHV1-1",
        ),
    ]
    return pd.DataFrame(rows)


def _base_tree(**overrides) -> AntibodyTree:
    t = AntibodyTree(
        family_id="d1::f1",
        root_id="root",
        node_ids=["root", "i", "l1"],
        edges=[("root", "i"), ("i", "l1")],
        branch_lengths={("root", "i"): 0.1, ("i", "l1"): 0.2},
        true_sequences={
            "root": "ATGCAT",
            "i": "ATGTAT",
            "l1": "ATCTAT",
        },
        is_leaf={"root": False, "i": False, "l1": True},
    )
    for k, v in overrides.items():
        setattr(t, k, v)
    return t


def test_validate_tree_valid():
    ok, reason = validate_tree(_base_tree())
    assert ok, reason


def test_validate_tree_multi_root():
    t = _base_tree(
        node_ids=["root", "i", "l1", "orphan"],
        true_sequences={
            "root": "ATGCAT",
            "i": "ATGTAT",
            "l1": "ATCTAT",
            "orphan": "AAAAAA",
        },
        is_leaf={"root": False, "i": False, "l1": True, "orphan": True},
    )
    ok, reason = validate_tree(t)
    assert not ok
    assert reason == "multi_root"


def test_validate_tree_disconnected():
    # Unique topo-root plus unreachable side component (b<->c).
    t = AntibodyTree(
        family_id="x",
        root_id="root",
        node_ids=["root", "a", "b", "c"],
        edges=[("root", "a"), ("b", "c"), ("c", "b")],
        branch_lengths={("root", "a"): 0.1, ("b", "c"): 0.1, ("c", "b"): 0.1},
        true_sequences={
            "root": "AAA",
            "a": "AAT",
            "b": "ATA",
            "c": "TTA",
        },
        is_leaf={"root": False, "a": True, "b": False, "c": False},
    )
    ok, reason = validate_tree(t)
    assert not ok
    assert reason == "disconnected"


def test_validate_tree_cycle():
    t = AntibodyTree(
        family_id="x",
        root_id="root",
        node_ids=["root", "a", "b"],
        edges=[("root", "a"), ("a", "b"), ("b", "a")],
        branch_lengths={("root", "a"): 0.1, ("a", "b"): 0.1, ("b", "a"): 0.1},
        true_sequences={"root": "AAA", "a": "AAT", "b": "ATA"},
        is_leaf={"root": False, "a": False, "b": False},
    )
    ok, reason = validate_tree(t)
    assert not ok
    assert reason == "cycle"


def test_validate_tree_missing_seq():
    t = _base_tree(true_sequences={"root": "ATGCAT", "i": "ATGTAT", "l1": ""})
    ok, reason = validate_tree(t)
    assert not ok
    assert reason == "missing_seq"


def test_validate_tree_invalid_bl():
    t = _base_tree(branch_lengths={("root", "i"): 0.0, ("i", "l1"): 0.2})
    ok, reason = validate_tree(t)
    assert not ok
    assert reason == "nonpos_bl"
    t2 = _base_tree(branch_lengths={("root", "i"): float("nan"), ("i", "l1"): 0.2})
    ok2, reason2 = validate_tree(t2)
    assert not ok2
    assert reason2 in {"nonpos_bl", "invalid_bl"}


def test_reconstruct_prefers_naive_root():
    trees, attrition, reasons, stages = reconstruct_trees_from_pcp_df(
        _toy_pcp_df(),
        min_leaves=4,
        min_depth=2,
        require_single_root=True,
        return_attrition=True,
    )
    assert len(trees) == 1
    assert trees[0].root_id == "root"
    assert trees[0].metadata["root_source"] == "naive"
    ok, reason = validate_tree(trees[0])
    assert ok, reason
    assert stages["final_eligible"] == 1


def test_reconstruct_excludes_multi_root_without_naive():
    rows = [
        dict(
            sample_id="d1",
            family="forest",
            parent_name="r1",
            child_name="c1",
            parent_heavy="ATGCATGCATGC",
            child_heavy="ATGTATGCATGC",
            branch_length=0.05,
            depth=1,
            parent_is_naive=False,
            child_is_leaf=True,
        ),
        dict(
            sample_id="d1",
            family="forest",
            parent_name="r2",
            child_name="c2",
            parent_heavy="ATGCATGCATGC",
            child_heavy="ATGCATGCATGT",
            branch_length=0.05,
            depth=1,
            parent_is_naive=False,
            child_is_leaf=True,
        ),
    ]
    trees, attrition, reasons, stages = reconstruct_trees_from_pcp_df(
        pd.DataFrame(rows),
        min_leaves=1,
        min_depth=1,
        return_attrition=True,
    )
    assert trees == []
    assert reasons.get("multi_root_no_usable_naive", 0) == 1


def test_rollout_api_hides_nonroot_truth_and_is_recursive():
    trees = reconstruct_trees_from_pcp_df(_toy_pcp_df(), min_leaves=4, min_depth=2)
    tree = trees[0]
    model = PoissonSiteMutationNullModel()
    root_aa = "AAAA"
    gen = rollout_tree(
        model,
        root_aa,
        tree.edges,
        tree.branch_lengths,
        seed=1,
        root_id=tree.root_id,
        family_id=tree.family_id,
    )
    assert gen[tree.root_id] == "AAAA"
    assert "i" in gen
    # Identity null keeps parent — verifies API path; no true_sequences involved
    idm = IdentityNullModel()
    gen2 = rollout_tree(
        idm,
        root_aa,
        tree.edges,
        tree.branch_lengths,
        seed=1,
        root_id=tree.root_id,
        family_id=tree.family_id,
    )
    assert all(v == "AAAA" for v in gen2.values())
