"""TreeSBM dataset: loads preprocessed trees and returns (root_seq, TreeState) pairs."""
import pickle
from pathlib import Path
from typing import Optional, Tuple

import torch
from torch.utils.data import Dataset

from src.tree_state import TreeState


class EvolutionaryTreeDataset(Dataset):
    """
    Dataset of evolutionary trees with root sequences.
    Each example is a (root_sequence, TreeState) pair from preprocessed .pkl files.
    """

    def __init__(
        self,
        processed_dir: str = "data/processed",
        split: str = "train",
        held_out_pattern: Optional[str] = None,
        test_size: float = 0.1,
    ):
        """
        Args:
            processed_dir: directory containing *_tree_data.pkl files from Phase 0
            split: "train" or "test"
            held_out_pattern: if set, trees matching this pattern go to test set
                e.g., "kinase" → all trees with "kinase" in name are test
            test_size: if held_out_pattern is None, use random split with this fraction
        """
        self.processed_dir = Path(processed_dir)
        self.split = split
        self.held_out_pattern = held_out_pattern
        self.test_size = test_size

        # Find all pickle files
        pkl_files = sorted(self.processed_dir.glob("*_tree_data.pkl"))
        if not pkl_files:
            raise ValueError(f"No *_tree_data.pkl files found in {processed_dir}")

        # Split into train/test
        if held_out_pattern:
            # Pattern-based split
            test_files = [
                f for f in pkl_files if held_out_pattern.lower() in f.stem.lower()
            ]
            train_files = [f for f in pkl_files if f not in test_files]
        else:
            # Random split
            n_test = int(len(pkl_files) * test_size)
            n_test = max(0, min(n_test, len(pkl_files) - 1))  # keep at least 1 for train
            import random

            random.seed(42)
            shuffled = random.sample(pkl_files, len(pkl_files))
            test_files = shuffled[:n_test]
            train_files = shuffled[n_test:]

        self.files = test_files if split == "test" else train_files

    def __len__(self) -> int:
        return len(self.files)

    def __getitem__(self, idx: int) -> Tuple[str, TreeState]:
        """
        Returns:
            (root_sequence, tree_state)
            - root_sequence: str, the root amino acid sequence
            - tree_state: TreeState object representing the full tree
        """
        pkl_path = self.files[idx]

        with open(pkl_path, "rb") as f:
            pkl_data = pickle.load(f)

        # Extract root sequence
        root_seq = pkl_data["root_seq"]

        # Create TreeState from pickle data
        tree_state = TreeState.from_newick_pkl(pkl_data)

        return root_seq, tree_state
