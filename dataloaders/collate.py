"""Custom collate function for batching variable-size trees."""
from typing import Any, Dict, List, Tuple

import numpy as np
import torch

from src.tree_state import TreeState


def collate_tree_batch(batch: List[Tuple[str, TreeState]]) -> Dict[str, Any]:
    """
    Collate function for batching trees with variable structure and sequence lengths.

    Input:
        batch: list of (root_seq: str, tree_state: TreeState) tuples

    Output:
        dict with keys:
            - root_seqs: (N, max_seq_len) tensor of root sequences, padded
            - node_seqs: (total_nodes, max_seq_len) tensor of all node sequences, padded
            - edge_index: (2, total_edges) tensor of graph connectivity
            - branch_lengths: (total_edges,) tensor of edge weights
            - graph_assignment: (total_nodes,) tensor mapping nodes to batch examples
            - active_leaves_mask: (total_nodes,) boolean tensor for leaf nodes
            - n_nodes_per_graph: list of node counts per tree (for unbatching later)
    """
    root_seqs, tree_states = zip(*batch)
    batch_size = len(batch)

    # Find max sequence length across all sequences in batch
    max_seq_len = max(
        max(len(root_seq), max([len(seq) for seq in ts.node_seqs.values()]))
        for root_seq, ts in zip(root_seqs, tree_states)
    )

    # Pad sequences and build node-level data
    padded_root_seqs = []
    padded_node_seqs_list = []
    node_counts = []
    graph_assignments = []
    active_leaves_masks = []

    node_offset = 0
    all_edges = []
    all_branch_lengths = []

    for graph_idx, (root_seq, tree_state) in enumerate(zip(root_seqs, tree_states)):
        # Pad root sequence
        padded_root = _pad_sequence(root_seq, max_seq_len)
        padded_root_seqs.append(padded_root)

        # Pad all node sequences for this tree
        padded_nodes = []
        for node_id in tree_state.node_ids:
            seq = tree_state.node_seqs[node_id]
            padded_seq = _pad_sequence(seq, max_seq_len)
            padded_nodes.append(padded_seq)

        padded_node_seqs_list.extend(padded_nodes)

        # Track graph assignment (which batch example each node belongs to)
        n_nodes = len(tree_state.node_ids)
        graph_assignments.extend([graph_idx] * n_nodes)
        node_counts.append(n_nodes)

        # Track active leaves (leaf nodes that are still growing)
        active_leaves = set(tree_state.active_leaves)
        for node_id in tree_state.node_ids:
            active_leaves_masks.append(node_id in active_leaves)

        # Build edge index with offset
        for parent_id, child_id in tree_state.edges:
            parent_idx = tree_state.node_ids.index(parent_id) + node_offset
            child_idx = tree_state.node_ids.index(child_id) + node_offset
            all_edges.append([parent_idx, child_idx])
            all_branch_lengths.append(
                tree_state.branch_lengths.get((parent_id, child_id), 0.0)
            )

        node_offset += n_nodes

    # Convert to tensors
    root_seqs_tensor = torch.tensor(
        padded_root_seqs, dtype=torch.long
    )  # (N, max_seq_len)
    node_seqs_tensor = torch.tensor(
        padded_node_seqs_list, dtype=torch.long
    )  # (total_nodes, max_seq_len)

    if all_edges:
        edge_index = torch.tensor(all_edges, dtype=torch.long).t()  # (2, total_edges)
        branch_lengths = torch.tensor(all_branch_lengths, dtype=torch.float32)
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)
        branch_lengths = torch.zeros(0, dtype=torch.float32)

    graph_assignment = torch.tensor(graph_assignments, dtype=torch.long)
    active_leaves_mask = torch.tensor(active_leaves_masks, dtype=torch.bool)

    return {
        "root_seqs": root_seqs_tensor,
        "node_seqs": node_seqs_tensor,
        "edge_index": edge_index,
        "branch_lengths": branch_lengths,
        "graph_assignment": graph_assignment,
        "active_leaves_mask": active_leaves_mask,
        "n_nodes_per_graph": node_counts,
    }


def _pad_sequence(seq: str, target_len: int, pad_token: str = "-") -> List[int]:
    """
    Convert amino acid sequence to padded integer indices.

    Amino acids are encoded as:
        A=0, C=1, D=2, ..., Y=24, -=25 (gap/pad)

    Args:
        seq: amino acid sequence string
        target_len: target length after padding
        pad_token: character to use for padding (default: '-' for gaps)

    Returns:
        list of integers of length target_len
    """
    # Standard amino acid alphabet + gap
    aa_to_idx = {
        "A": 0,
        "C": 1,
        "D": 2,
        "E": 3,
        "F": 4,
        "G": 5,
        "H": 6,
        "I": 7,
        "K": 8,
        "L": 9,
        "M": 10,
        "N": 11,
        "P": 12,
        "Q": 13,
        "R": 14,
        "S": 15,
        "T": 16,
        "V": 17,
        "W": 18,
        "Y": 19,
        "-": 20,  # gap
        "X": 20,  # unknown → treat as gap
    }
    pad_idx = 20  # index for padding

    # Convert to indices
    indices = [aa_to_idx.get(aa.upper(), pad_idx) for aa in seq]

    # Pad to target length
    if len(indices) < target_len:
        indices.extend([pad_idx] * (target_len - len(indices)))
    else:
        indices = indices[:target_len]

    return indices
