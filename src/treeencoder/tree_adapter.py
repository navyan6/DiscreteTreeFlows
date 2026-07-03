from collections.abc import Callable

import torch

from src.tree_state import TreeState
from src.treeencoder.encoder_input import TreeEncoderInput
from src.treeencoder.structural_features import compute_structural_features
from src.treeencoder.edges import build_edges
from src.treeencoder.laplacian import compute_laplacian_pe

def tree_state_to_encoder_input(
        tree: TreeState, 
        node_embeddings: torch.Tensor, 
        laplacian_dim: int, 
) -> TreeEncoderInput: 

    """
    Take a treestate and return a tensor representation of the treestate
    """

    num_nodes = tree.n_nodes()

    if node_embeddings.ndim != 2:
        raise ValueError("node embeddings should have 2 values")
    
    if node_embeddings.shape[0] != num_nodes:
        raise ValueError("there should be one row per node")
    
    node_to_index = {node_id:index for index, node_id in enumerate(tree.node_ids)}
    
    structural_features = compute_structural_features(
        tree=tree,
        node_to_index=node_to_index,
    )

    edge_idx, edge_types, edge_attr = build_edges(
        tree=tree,
        node_to_idx=node_to_index,
    )

    lap_pe = compute_laplacian_pe(
        tree=tree,
        node_to_index=node_to_index,
        num_eigenvectors=laplacian_dim,
    )

    return TreeEncoderInput(
        node_ids=tree.node_ids,
        x=node_embeddings,
        structural_features=structural_features,
        lap_pe=lap_pe,
        edge_index=edge_idx,
        edge_type=edge_types,
        edge_attr=edge_attr,
        root_index=node_to_index[tree.root_id],
    )

