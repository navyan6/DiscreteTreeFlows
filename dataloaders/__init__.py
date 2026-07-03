"""Data loading for TreeSBM."""
from .tree_dataset import EvolutionaryTreeDataset
from .collate import collate_tree_batch

__all__ = ["EvolutionaryTreeDataset", "collate_tree_batch"]
