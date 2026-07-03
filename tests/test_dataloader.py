"""Tests for DataLoader components."""
import pickle
import tempfile
from pathlib import Path

import torch

from dataloaders.tree_dataset import EvolutionaryTreeDataset
from dataloaders.collate import collate_tree_batch, _pad_sequence
from src.tree_state import TreeState


class TestPadSequence:
    """Test sequence padding."""

    def test_pad_short_sequence(self):
        """Pad short sequence."""
        indices = _pad_sequence("MK", 5)
        assert len(indices) == 5
        assert indices[0] == 10  # M
        assert indices[1] == 8  # K
        assert indices[2:] == [20, 20, 20]  # padding

    def test_no_padding_needed(self):
        """Exact length sequence."""
        indices = _pad_sequence("MKVL", 4)
        assert len(indices) == 4
        assert indices == [10, 8, 17, 9]  # M, K, V, L

    def test_truncate_long_sequence(self):
        """Truncate longer sequence."""
        indices = _pad_sequence("MKVLAA", 4)
        assert len(indices) == 4
        assert indices == [10, 8, 17, 9]  # M, K, V, L (truncated)

    def test_gap_character(self):
        """Handle gap character."""
        indices = _pad_sequence("M-KV", 4)
        assert indices[1] == 20  # gap → index 20


class TestDataset:
    """Test EvolutionaryTreeDataset."""

    @staticmethod
    def _create_test_pkl(tree_name: str, pkl_dir: Path):
        """Create a test tree pickle file."""
        pkl_data = {
            "name": tree_name,
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
        pkl_path = pkl_dir / f"{tree_name}_tree_data.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(pkl_data, f)

    def test_dataset_length(self):
        """Dataset reports correct length."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            for i in range(5):
                self._create_test_pkl(f"tree_{i}", tmpdir)

            dataset = EvolutionaryTreeDataset(tmpdir)
            assert len(dataset) == 5

    def test_dataset_getitem(self):
        """Get item returns (root_seq, TreeState)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._create_test_pkl("test_tree", tmpdir)

            dataset = EvolutionaryTreeDataset(tmpdir)
            root_seq, tree_state = dataset[0]

            assert isinstance(root_seq, str)
            assert isinstance(tree_state, TreeState)
            assert root_seq == "MKVL"
            assert tree_state.n_nodes() == 3

    def test_train_test_split(self):
        """Random train/test split works."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            for i in range(10):
                self._create_test_pkl(f"tree_{i}", tmpdir)

            train_ds = EvolutionaryTreeDataset(tmpdir, split="train", test_size=0.2)
            test_ds = EvolutionaryTreeDataset(tmpdir, split="test", test_size=0.2)

            assert len(train_ds) == 8
            assert len(test_ds) == 2
            assert len(train_ds) + len(test_ds) == 10

    def test_pattern_based_split(self):
        """Pattern-based train/test split."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            for i in range(5):
                self._create_test_pkl(f"kinase_{i}", tmpdir)
            for i in range(5):
                self._create_test_pkl(f"phosphatase_{i}", tmpdir)

            train_ds = EvolutionaryTreeDataset(
                tmpdir, split="train", held_out_pattern="kinase"
            )
            test_ds = EvolutionaryTreeDataset(
                tmpdir, split="test", held_out_pattern="kinase"
            )

            assert len(train_ds) == 5  # phosphatases
            assert len(test_ds) == 5  # kinases


class TestCollate:
    """Test collate_tree_batch."""

    @staticmethod
    def _create_sample(root_seq: str, n_children: int = 2):
        """Create a sample (root_seq, TreeState) tuple."""
        tree = TreeState.root_only(root_seq)
        child_seqs = [root_seq.replace("M", "A") for _ in range(n_children)]
        tree = tree.branch_node("root", child_seqs)
        return root_seq, tree

    def test_collate_single_tree(self):
        """Collate single tree."""
        batch = [self._create_sample("MK")]
        output = collate_tree_batch(batch)

        assert output["root_seqs"].shape == (1, 2)
        assert output["node_seqs"].shape[0] == 3  # root + 2 children
        assert output["edge_index"].shape[1] == 2  # 2 edges
        assert output["graph_assignment"].shape[0] == 3

    def test_collate_variable_lengths(self):
        """Collate trees with different sequence lengths."""
        batch = [
            self._create_sample("MK"),
            self._create_sample("MKVL"),
            self._create_sample("MKVLAA"),
        ]
        output = collate_tree_batch(batch)

        # All sequences padded to max length
        assert output["root_seqs"].shape == (3, 6)
        assert output["node_seqs"].shape[1] == 6

    def test_collate_batch_size(self):
        """Collate multiple trees."""
        batch = [self._create_sample("MK") for _ in range(4)]
        output = collate_tree_batch(batch)

        assert output["root_seqs"].shape[0] == 4
        assert output["graph_assignment"].nunique() == 4  # 4 different graphs

    def test_collate_edge_offset(self):
        """Edge indices are correctly offset per graph."""
        batch = [self._create_sample("MK"), self._create_sample("KV")]
        output = collate_tree_batch(batch)

        edge_index = output["edge_index"]
        # First tree edges: (0, 1), (0, 2)
        # Second tree edges: (3, 4), (3, 5) (offset by 3)
        assert edge_index.shape[1] == 4
        assert edge_index[0, 0].item() == 0  # first edge parent is node 0
        assert edge_index[0, 2].item() == 3  # third edge parent is node 3 (offset)


if __name__ == "__main__":
    import pytest

    pytest.main([__file__, "-v"])
