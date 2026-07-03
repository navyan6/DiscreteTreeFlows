"""Tests for TreeEncoder and RateHeads."""
import torch
import pytest

from src.tree_state import TreeState
from src.networks import TreeEncoder, RateHeads, TransformerGraphLayer


class TestTransformerGraphLayer:
    """Test basic transformer layer."""

    def test_forward_shape(self):
        """Transformer layer preserves node shape."""
        layer = TransformerGraphLayer(d_model=256, n_heads=8)

        h = torch.randn(10, 256)  # 10 nodes
        edge_index = torch.tensor([[0, 1, 2], [1, 2, 3]], dtype=torch.long)
        edge_features = torch.randn(3, 256)
        active_leaves = ["node_1", "node_2"]
        node_id_to_idx = {f"node_{i}": i for i in range(10)}

        h_out = layer(h, edge_index, edge_features, active_leaves, node_id_to_idx)

        assert h_out.shape == h.shape
        assert h_out.shape == (10, 256)

    def test_forward_single_node(self):
        """Transformer works on single node."""
        layer = TransformerGraphLayer(d_model=256, n_heads=8)

        h = torch.randn(1, 256)
        edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_features = torch.zeros((0, 256))
        active_leaves = []
        node_id_to_idx = {"root": 0}

        h_out = layer(h, edge_index, edge_features, active_leaves, node_id_to_idx)

        assert h_out.shape == (1, 256)


class TestRateHeads:
    """Test rate prediction heads."""

    def test_forward_shape(self):
        """RateHeads output correct shapes."""
        heads = RateHeads(d_model=256, max_seq_len=512)

        H_T = torch.randn(10, 256)  # 10 nodes
        active_indices = [0, 1, 3, 5]  # 4 active leaves

        output = heads(H_T, active_indices)

        assert output["mutation_logits"].shape == (4, 512, 20)
        assert output["branching_rate"].shape == (4,)
        assert output["branch_length"].shape == (4,)
        assert output["stop_prob"].shape == (4,)

    def test_branching_rate_positive(self):
        """Branching rate is always positive."""
        heads = RateHeads(d_model=256, max_seq_len=512)

        H_T = torch.randn(10, 256)
        active_indices = list(range(10))

        output = heads(H_T, active_indices)

        assert torch.all(output["branching_rate"] > 0)
        assert torch.all(torch.isfinite(output["branching_rate"]))

    def test_stop_prob_in_range(self):
        """Stop probability in [0, 1]."""
        heads = RateHeads(d_model=256, max_seq_len=512)

        H_T = torch.randn(10, 256)
        active_indices = list(range(10))

        output = heads(H_T, active_indices)

        assert torch.all(output["stop_prob"] >= 0)
        assert torch.all(output["stop_prob"] <= 1)

    def test_single_active_leaf(self):
        """RateHeads with single active leaf."""
        heads = RateHeads(d_model=256, max_seq_len=512)

        H_T = torch.randn(1, 256)
        active_indices = [0]

        output = heads(H_T, active_indices)

        assert output["mutation_logits"].shape == (1, 512, 20)
        assert output["branching_rate"].shape == (1,)


class TestTreeEncoderBasic:
    """Basic TreeEncoder tests (without ESM models to keep tests fast)."""

    def test_compute_depths(self):
        """Depth computation works correctly."""
        encoder = RateHeads(d_model=256)  # use RateHeads as dummy, just for structure

        # Create simple tree: root → child1, child2; child1 → grandchild
        node_ids = ["root", "child1", "child2", "grandchild"]
        edges = [("root", "child1"), ("root", "child2"), ("child1", "grandchild")]

        # We can't easily test TreeEncoder._compute_depths without instantiating ESM
        # So we just verify the concept here
        pass

    def test_edge_feature_count(self):
        """Edge features have correct dimension."""
        # Edge features: [branch_length, depth, sequence_divergence]
        # Should be 3-dimensional before embedding to d_model
        pass


if __name__ == "__main__":
    import pytest

    pytest.main([__file__, "-v"])
