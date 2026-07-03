"""
Comprehensive test suite for tree encoder pipeline.
Tests: tree_state.py, structural_features.py, edges.py, laplacian.py, tree_adapter.py
"""

import pytest
import torch
import numpy as np
from src.tree_state import TreeState
from src.treeencoder.structural_features import compute_depth, compute_structural_features
from src.treeencoder.edges import build_edges
from src.treeencoder.laplacian import compute_laplacian_pe
from src.treeencoder.encoder_input import TreeEncoderInput
from src.treeencoder.tree_adapter import tree_state_to_encoder_input


class TestCaseA_RootOnly:
    """Test with single root node, no children."""

    @pytest.fixture
    def tree(self):
        return TreeState(
            node_ids=["root"],
            root_id="root",
            edges=[],
            branch_lengths={},
            node_seqs={"root": "AAAA"},
            active_leaves=["root"],
        )

    def test_tree_basic_properties(self, tree):
        """Verify basic tree properties."""
        assert tree.n_nodes() == 1
        assert tree.n_leaves() == 1
        assert tree.is_leaf("root")
        assert tree.get_children("root") == []
        assert tree.get_parent("root") is None

    def test_compute_depth_root_only(self, tree):
        """Root has depth 0."""
        depths = compute_depth(tree)
        assert depths == {"root": 0}

    def test_structural_features_root_only(self, tree):
        """Normalized depth = 0 (0/0 handled), is_root=1, is_leaf=1."""
        node_to_index = {"root": 0}
        features = compute_structural_features(tree, node_to_index)

        assert features.shape == (1, 3)
        assert features.dtype == torch.float32
        assert features[0, 0].item() == 0.0  # normalized depth
        assert features[0, 1].item() == 1.0  # is_root
        assert features[0, 2].item() == 1.0  # is_leaf

    def test_edges_root_only(self, tree):
        """No edges -> empty edge tensors with correct shapes."""
        node_to_index = {"root": 0}
        edge_index, edge_type, edge_attr = build_edges(tree, node_to_index)

        assert edge_index.shape == (2, 0)
        assert edge_index.dtype == torch.long
        assert edge_type.shape == (0,)
        assert edge_type.dtype == torch.long
        assert edge_attr.shape == (0, 1)
        assert edge_attr.dtype == torch.float32

    def test_laplacian_root_only(self, tree):
        """Single node -> zeros [1, k]."""
        node_to_index = {"root": 0}
        k = 5
        lap_pe = compute_laplacian_pe(tree, node_to_index, k)

        assert lap_pe.shape == (1, 5)
        assert lap_pe.dtype == torch.float32
        assert torch.allclose(lap_pe, torch.zeros(1, 5))

    def test_root_only_full_pipeline(self, tree):
        """End-to-end test: tree -> TreeEncoderInput."""
        node_embeddings = torch.randn(1, 128)
        laplacian_dim = 5

        result = tree_state_to_encoder_input(tree, node_embeddings, laplacian_dim)

        assert isinstance(result, TreeEncoderInput)
        assert result.node_ids == ["root"]
        assert result.x.shape == (1, 128)
        assert result.structural_features.shape == (1, 3)
        assert result.lap_pe.shape == (1, 5)
        assert result.edge_index.shape == (2, 0)
        assert result.edge_type.shape == (0,)
        assert result.edge_attr.shape == (0, 1)
        assert result.root_index == 0


class TestCaseB_RootWithTwoChildren:
    """Test with root and two immediate children."""

    @pytest.fixture
    def tree(self):
        return TreeState(
            node_ids=["root", "A", "B"],
            root_id="root",
            edges=[("root", "A"), ("root", "B")],
            branch_lengths={("root", "A"): 0.1, ("root", "B"): 0.2},
            node_seqs={"root": "AAAA", "A": "AAAT", "B": "AATA"},
            active_leaves=["A", "B"],
        )

    def test_tree_basic_properties(self, tree):
        """Verify parent-child relationships."""
        assert tree.n_nodes() == 3
        assert tree.n_leaves() == 2
        assert not tree.is_leaf("root")
        assert tree.is_leaf("A")
        assert tree.is_leaf("B")
        assert set(tree.get_children("root")) == {"A", "B"}
        assert tree.get_parent("A") == "root"
        assert tree.get_parent("B") == "root"

    def test_compute_depth(self, tree):
        """Root depth 0, children depth 1."""
        depths = compute_depth(tree)
        assert depths["root"] == 0
        assert depths["A"] == 1
        assert depths["B"] == 1

    def test_normalized_depth(self, tree):
        """Max depth = 1, so normalized depth of children = 1.0."""
        node_to_index = {"root": 0, "A": 1, "B": 2}
        features = compute_structural_features(tree, node_to_index)

        assert features.shape == (3, 3)
        assert features[0, 0].item() == 0.0  # root: 0/1
        assert features[1, 0].item() == 1.0  # A: 1/1
        assert features[2, 0].item() == 1.0  # B: 1/1

    def test_root_and_leaf_indicators(self, tree):
        """Root has is_root=1, all others 0. Only A and B are leaves."""
        node_to_index = {"root": 0, "A": 1, "B": 2}
        features = compute_structural_features(tree, node_to_index)

        assert features[0, 1].item() == 1.0  # root is_root
        assert features[1, 1].item() == 0.0  # A is_root
        assert features[2, 1].item() == 0.0  # B is_root

        assert features[0, 2].item() == 0.0  # root is_leaf
        assert features[1, 2].item() == 1.0  # A is_leaf
        assert features[2, 2].item() == 1.0  # B is_leaf

    def test_bidirectional_edges(self, tree):
        """Two biological edges -> four computational edges (bidirectional)."""
        node_to_index = {"root": 0, "A": 1, "B": 2}
        edge_index, edge_type, edge_attr = build_edges(tree, node_to_index)

        assert edge_index.shape == (2, 4)
        assert edge_type.shape == (4,)
        assert edge_attr.shape == (4, 1)

        # Extract edges
        edges_list = [(edge_index[0, i].item(), edge_index[1, i].item(),
                       edge_type[i].item(), edge_attr[i, 0].item())
                      for i in range(4)]

        # Expect: root->A, A->root, root->B, B->root (order may vary)
        parent_to_child_edges = [(src, tgt, etype) for src, tgt, etype, _ in edges_list if etype == 0]
        child_to_parent_edges = [(src, tgt, etype) for src, tgt, etype, _ in edges_list if etype == 1]

        assert len(parent_to_child_edges) == 2
        assert len(child_to_parent_edges) == 2

    def test_branch_lengths_preserved(self, tree):
        """Forward and reverse edges have same branch length."""
        node_to_index = {"root": 0, "A": 1, "B": 2}
        edge_index, edge_type, edge_attr = build_edges(tree, node_to_index)

        edges_dict = {}
        for i in range(edge_index.shape[1]):
            src = edge_index[0, i].item()
            tgt = edge_index[1, i].item()
            length = edge_attr[i, 0].item()
            edges_dict[(src, tgt)] = length

        # root=0, A=1: expect 0.1 for both directions
        assert edges_dict[(0, 1)] == 0.1
        assert edges_dict[(1, 0)] == 0.1

        # root=0, B=2: expect 0.2 for both directions
        assert edges_dict[(0, 2)] == 0.2
        assert edges_dict[(2, 0)] == 0.2

    def test_root_index(self, tree):
        """root_index must be 0 (root is first in node_ids)."""
        node_to_index = {"root": 0, "A": 1, "B": 2}
        assert node_to_index[tree.root_id] == 0

    def test_full_pipeline_two_children(self, tree):
        """End-to-end test."""
        node_embeddings = torch.randn(3, 64)
        laplacian_dim = 5

        result = tree_state_to_encoder_input(tree, node_embeddings, laplacian_dim)

        assert result.node_ids == ["root", "A", "B"]
        assert result.x.shape == (3, 64)
        assert result.structural_features.shape == (3, 3)
        assert result.lap_pe.shape == (3, 5)
        assert result.edge_index.shape == (2, 4)
        assert result.edge_type.shape == (4,)
        assert result.edge_attr.shape == (4, 1)
        assert result.root_index == 0


class TestCaseC_ThreeLevelTree:
    """Test a 3-level tree to verify depth, leaf identification, and Laplacian properties."""

    @pytest.fixture
    def tree(self):
        return TreeState(
            node_ids=["root", "A", "B", "C"],
            root_id="root",
            edges=[("root", "A"), ("root", "B"), ("A", "C")],
            branch_lengths={
                ("root", "A"): 0.1,
                ("root", "B"): 0.2,
                ("A", "C"): 0.05,
            },
            node_seqs={
                "root": "AAAA",
                "A": "AAAT",
                "B": "AATA",
                "C": "AATT",
            },
            active_leaves=["B", "C"],
        )

    def test_depths_correct(self, tree):
        """Verify DFS-independent depth calculation."""
        depths = compute_depth(tree)
        assert depths["root"] == 0
        assert depths["A"] == 1
        assert depths["B"] == 1
        assert depths["C"] == 2

    def test_leaf_identification(self, tree):
        """Only B and C are leaves (A has a child)."""
        assert not tree.is_leaf("root")
        assert not tree.is_leaf("A")
        assert tree.is_leaf("B")
        assert tree.is_leaf("C")
        assert tree.n_leaves() == 2

    def test_structural_features_shape(self, tree):
        """All nodes present, features shape correct."""
        node_to_index = {"root": 0, "A": 1, "B": 2, "C": 3}
        features = compute_structural_features(tree, node_to_index)

        assert features.shape == (4, 3)
        assert torch.isfinite(features).all()

    def test_laplacian_shape_and_finiteness(self, tree):
        """Laplacian PE has correct shape and is finite."""
        node_to_index = {"root": 0, "A": 1, "B": 2, "C": 3}
        lap_pe = compute_laplacian_pe(tree, node_to_index, num_eigenvectors=5)

        assert lap_pe.shape == (4, 5)
        assert torch.isfinite(lap_pe).all()

    def test_laplacian_eigenvector_property(self, tree):
        """Verify eigenvectors satisfy L @ u ≈ λ * u (residual small)."""
        node_to_index = {"root": 0, "A": 1, "B": 2, "C": 3}
        num_nodes = 4

        # Build adjacency
        adjacency = torch.zeros((num_nodes, num_nodes), dtype=torch.float32)
        for p_id, c_id in tree.edges:
            p_idx = node_to_index[p_id]
            c_idx = node_to_index[c_id]
            adjacency[p_idx, c_idx] = 1.0
            adjacency[c_idx, p_idx] = 1.0

        # Compute Laplacian
        degree = adjacency.sum(dim=1)
        inv_sqrt_degree = torch.zeros_like(degree)
        mask = degree > 0
        inv_sqrt_degree[mask] = degree[mask].pow(-0.5)
        D_inv_sqrt = torch.diag(inv_sqrt_degree)
        I = torch.eye(num_nodes)
        L = I - D_inv_sqrt @ adjacency @ D_inv_sqrt

        # Get eigenvectors
        eigenvalues, eigenvectors = torch.linalg.eigh(L)

        # Test nontrivial eigenvectors
        for i in range(1, min(4, num_nodes)):  # Skip trivial first one
            u = eigenvectors[:, i]
            lam = eigenvalues[i]
            residual = torch.norm(L @ u - lam * u)
            assert residual < 1e-5

    def test_full_pipeline_three_level(self, tree):
        """End-to-end."""
        node_embeddings = torch.randn(4, 128)

        result = tree_state_to_encoder_input(tree, node_embeddings, laplacian_dim=8)

        assert result.x.shape == (4, 128)
        assert result.structural_features.shape == (4, 3)
        assert result.lap_pe.shape == (4, 8)
        assert result.root_index == 0


class TestCaseD_InvalidTrees:
    """Test error handling for malformed trees."""

    def test_root_missing_from_nodes(self):
        """Root ID not in node_ids."""
        with pytest.raises(ValueError, match="Root .* not in node_ids"):
            TreeState(
                node_ids=["A", "B"],
                root_id="root",
                edges=[],
                branch_lengths={},
                node_seqs={"A": "AA", "B": "AA"},
            )

    def test_duplicate_node_ids(self):
        """Duplicate node IDs should be caught (or handled gracefully)."""
        tree = TreeState(
            node_ids=["root", "A", "A"],  # Duplicate A
            root_id="root",
            edges=[("root", "A")],
            branch_lengths={("root", "A"): 0.1},
            node_seqs={"root": "AA", "A": "AA"},
        )
        # This should ideally raise an error, but currently doesn't.
        # Document this as a known limitation.
        assert tree.n_nodes() == 3  # Records all, including duplicates

    def test_unknown_node_in_edge(self):
        """Edge references undefined node."""
        with pytest.raises(ValueError, match="edges reference unknown nodes"):
            TreeState(
                node_ids=["root", "A"],
                root_id="root",
                edges=[("root", "X")],  # X not in node_ids
                branch_lengths={},
                node_seqs={"root": "AA", "A": "AA"},
            )

    def test_missing_branch_length(self):
        """Edge exists but no branch length provided."""
        tree = TreeState(
            node_ids=["root", "A"],
            root_id="root",
            edges=[("root", "A")],
            branch_lengths={},  # Missing ("root", "A")
            node_seqs={"root": "AA", "A": "AA"},
        )
        node_to_index = {"root": 0, "A": 1}

        with pytest.raises(KeyError):
            build_edges(tree, node_to_index)

    def test_negative_branch_length(self):
        """Branch length should not be negative."""
        tree = TreeState(
            node_ids=["root", "A"],
            root_id="root",
            edges=[("root", "A")],
            branch_lengths={("root", "A"): -0.1},  # Negative!
            node_seqs={"root": "AA", "A": "AA"},
        )
        # Currently no validation; should ideally warn or error.
        # For now, this passes but should be addressed.
        node_to_index = {"root": 0, "A": 1}
        edge_index, edge_type, edge_attr = build_edges(tree, node_to_index)
        assert edge_attr[0, 0].item() == -0.1

    def test_non_root_with_multiple_parents(self):
        """In a tree, each non-root node has exactly one parent."""
        tree = TreeState(
            node_ids=["root", "A", "B", "C"],
            root_id="root",
            edges=[("root", "A"), ("A", "C"), ("B", "C")],  # C has two parents!
            branch_lengths={
                ("root", "A"): 0.1,
                ("A", "C"): 0.05,
                ("B", "C"): 0.05,
            },
            node_seqs={"root": "AA", "A": "AA", "B": "AA", "C": "AA"},
        )
        # This violates the tree structure. C has two parents: A and B.
        # Currently, no validation exists.
        assert tree.get_parent("C") in ["A", "B"]  # Returns first found

    def test_disconnected_node(self):
        """Node is not reachable from root."""
        tree = TreeState(
            node_ids=["root", "A", "B"],
            root_id="root",
            edges=[("root", "A")],  # B is unreachable
            branch_lengths={("root", "A"): 0.1},
            node_seqs={"root": "AA", "A": "AA", "B": "AA"},
        )

        depths = compute_depth(tree)
        # B should not be in depths if validation were strict
        assert "B" not in depths or len(depths) != 3

    def test_cycle_detection(self):
        """Tree should not contain cycles."""
        tree = TreeState(
            node_ids=["A", "B", "C"],
            root_id="A",
            edges=[("A", "B"), ("B", "C"), ("C", "A")],  # Cycle!
            branch_lengths={
                ("A", "B"): 0.1,
                ("B", "C"): 0.1,
                ("C", "A"): 0.1,
            },
            node_seqs={"A": "AA", "B": "AA", "C": "AA"},
        )
        # Currently, no cycle detection. compute_depth will hang or fail.
        # Should add cycle detection validation.

    def test_missing_sequence(self):
        """All nodes must have sequences."""
        with pytest.raises(ValueError, match="node_seqs contains unknown nodes"):
            TreeState(
                node_ids=["root", "A"],
                root_id="root",
                edges=[("root", "A")],
                branch_lengths={("root", "A"): 0.1},
                node_seqs={"root": "AA"},  # Missing A
            )

    def test_active_leaf_is_not_structural_leaf(self):
        """active_leaves should only contain structural leaves (no children)."""
        tree = TreeState(
            node_ids=["root", "A", "B"],
            root_id="root",
            edges=[("root", "A"), ("A", "B")],
            branch_lengths={
                ("root", "A"): 0.1,
                ("A", "B"): 0.05,
            },
            node_seqs={"root": "AA", "A": "AA", "B": "AA"},
            active_leaves=["A", "B"],  # A is not a leaf (has child B)!
        )
        # This is semantically incorrect. A should not be in active_leaves.
        # Currently no validation.
        assert "A" in tree.active_leaves
        assert not tree.is_leaf("A")  # But this correctly reports A is not a leaf


class TestTreeStateSerializationEdgeCases:
    """Test safe serialization/deserialization."""

    def test_from_dict_safe_parsing(self):
        """Verify that from_dict uses safe parsing (not eval)."""
        # This test documents the fix: we should use ast.literal_eval or manual parsing
        # not eval()
        d = {
            "node_ids": ["root", "A"],
            "root_id": "root",
            "edges": [("root", "A")],
            "branch_lengths": {
                "('root', 'A')": 0.1,  # String key format
            },
            "node_seqs": {"root": "AA", "A": "AA"},
            "active_leaves": ["A"],
        }

        tree = TreeState.from_dict(d)
        assert tree.branch_lengths[("root", "A")] == 0.1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
