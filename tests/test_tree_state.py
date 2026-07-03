"""Tests for TreeState dataclass."""
import pytest
from src.tree_state import TreeState


class TestTreeStateBasics:
    """Test basic TreeState creation and properties."""

    def test_root_only(self):
        """Create single-node tree."""
        x0 = "MKVL"
        t = TreeState.root_only(x0)
        assert t.n_nodes() == 1
        assert t.n_leaves() == 1
        assert t.root_id == "root"
        assert t.node_seqs["root"] == x0
        assert t.active_leaves == ["root"]

    def test_root_only_long_sequence(self):
        """Create root with realistic protein sequence."""
        x0 = "M" + "K" * 1000
        t = TreeState.root_only(x0)
        assert t.n_nodes() == 1
        assert len(t.node_seqs["root"]) == 1001
        assert t.is_leaf("root")

    def test_get_children_leaf(self):
        """Leaf node has no children."""
        t = TreeState.root_only("MK")
        assert t.get_children("root") == []

    def test_get_parent_root(self):
        """Root has no parent."""
        t = TreeState.root_only("MK")
        assert t.get_parent("root") is None


class TestMutations:
    """Test mutation operations."""

    def test_apply_mutation(self):
        """Apply single amino acid mutation."""
        t = TreeState.root_only("MKVL")
        t2 = t.apply_mutation("root", 1, "A")
        assert t.node_seqs["root"] == "MKVL"  # original unchanged
        assert t2.node_seqs["root"] == "MAVL"
        assert t2.n_nodes() == 1

    def test_apply_mutation_multiple(self):
        """Apply multiple mutations sequentially."""
        t = TreeState.root_only("MKVL")
        t = t.apply_mutation("root", 0, "A")
        t = t.apply_mutation("root", 3, "A")
        assert t.node_seqs["root"] == "AKVA"

    def test_apply_mutation_invalid_position(self):
        """Reject mutation at invalid position."""
        t = TreeState.root_only("MK")
        with pytest.raises(ValueError):
            t.apply_mutation("root", 10, "A")

    def test_apply_mutation_invalid_node(self):
        """Reject mutation on nonexistent node."""
        t = TreeState.root_only("MK")
        with pytest.raises(ValueError):
            t.apply_mutation("nonexistent", 0, "A")


class TestBranching:
    """Test branching operations."""

    def test_branch_node_single_child(self):
        """Create single child."""
        t = TreeState.root_only("MKVL")
        t2 = t.branch_node("root", ["MKVL"])
        assert t2.n_nodes() == 2
        assert t2.n_leaves() == 1
        assert "root_child_0" in t2.node_ids
        assert ("root", "root_child_0") in t2.edges
        assert t2.node_seqs["root_child_0"] == "MKVL"
        assert t2.active_leaves == ["root_child_0"]

    def test_branch_node_multiple_children(self):
        """Create multiple children."""
        t = TreeState.root_only("MK")
        t2 = t.branch_node("root", ["MA", "ML", "MV"])
        assert t2.n_nodes() == 4
        assert t2.n_leaves() == 3
        assert t2.node_seqs["root_child_0"] == "MA"
        assert t2.node_seqs["root_child_1"] == "ML"
        assert t2.node_seqs["root_child_2"] == "MV"
        assert len(t2.active_leaves) == 3

    def test_branch_node_sequence_length_mismatch(self):
        """Reject children with wrong sequence length."""
        t = TreeState.root_only("MKVL")
        with pytest.raises(ValueError):
            t.branch_node("root", ["MK"])  # too short

    def test_branch_node_no_children(self):
        """Reject branching with no children."""
        t = TreeState.root_only("MK")
        with pytest.raises(ValueError):
            t.branch_node("root", [])


class TestBranchExtension:
    """Test branch length operations."""

    def test_extend_branch(self):
        """Extend branch leading to child."""
        t = TreeState.root_only("MK")
        t = t.branch_node("root", ["MA"])
        t2 = t.extend_branch("root_child_0", 0.5)
        assert t2.branch_lengths[("root", "root_child_0")] == 0.5

    def test_extend_branch_multiple(self):
        """Extend branch multiple times."""
        t = TreeState.root_only("MK")
        t = t.branch_node("root", ["MA"])
        t = t.extend_branch("root_child_0", 0.3)
        t = t.extend_branch("root_child_0", 0.2)
        assert t.branch_lengths[("root", "root_child_0")] == 0.5

    def test_extend_branch_root_error(self):
        """Cannot extend branch to root."""
        t = TreeState.root_only("MK")
        with pytest.raises(ValueError):
            t.extend_branch("root", 0.5)


class TestTermination:
    """Test leaf termination."""

    def test_terminate_leaf(self):
        """Terminate a leaf node."""
        t = TreeState.root_only("MK")
        t = t.branch_node("root", ["MA", "ML"])
        assert len(t.active_leaves) == 2
        t2 = t.terminate_leaf("root_child_0")
        assert "root_child_0" not in t2.active_leaves
        assert "root_child_1" in t2.active_leaves

    def test_terminate_non_active_error(self):
        """Cannot terminate non-active node."""
        t = TreeState.root_only("MK")
        with pytest.raises(ValueError):
            t.terminate_leaf("root")


class TestSerialization:
    """Test to_dict and from_dict."""

    def test_roundtrip_root_only(self):
        """Serialize and deserialize root-only tree."""
        x0 = "MKVL"
        t1 = TreeState.root_only(x0)
        d = t1.to_dict()
        t2 = TreeState.from_dict(d)
        assert t2.n_nodes() == t1.n_nodes()
        assert t2.root_id == t1.root_id
        assert t2.node_seqs == t1.node_seqs

    def test_roundtrip_with_branches(self):
        """Serialize and deserialize tree with branches."""
        t1 = TreeState.root_only("MK")
        t1 = t1.branch_node("root", ["MA", "ML"])
        t1 = t1.extend_branch("root_child_0", 0.5)
        d = t1.to_dict()
        t2 = TreeState.from_dict(d)
        assert t2.n_nodes() == 3
        assert t2.n_leaves() == 2
        assert t2.branch_lengths == t1.branch_lengths


class TestFromNewickPkl:
    """Test loading from preprocessed tree pickle."""

    def test_from_pkl_simple(self):
        """Load from Phase 0 pickle format."""
        pkl_data = {
            "name": "test_tree",
            "root_id": "root",
            "root_seq": "MKVL",
            "node_seqs": {
                "root": "MKVL",
                "leaf1": "MKVL",
                "leaf2": "MAVL",
            },
            "edges": [("root", "leaf1"), ("root", "leaf2")],
            "branch_lengths": {("root", "leaf1"): 0.1, ("root", "leaf2"): 0.2},
            "n_leaves": 2,
            "n_nodes": 3,
        }
        t = TreeState.from_newick_pkl(pkl_data)
        assert t.n_nodes() == 3
        assert t.n_leaves() == 2
        assert t.root_id == "root"
        assert len(t.active_leaves) == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
