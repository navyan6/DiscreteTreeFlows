#preserve edges in the graph
#reverse edges if necessary 

import torch

from src.tree_state import TreeState

parent_to_child = 0
child_to_parent = 1

#returns 3 tensors, of edge index, edge type, edge attribution
#edge_attr = branch length (sub/site), and we save source and target nodes
#add forward and reverse edges 

def build_edges(
    tree: TreeState,
    node_to_idx: dict[str, int],
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:

    sources: list[int] = []

    targets: list[int] = []

    edge_types: list[int] = []

    edge_attr: list[list[float]] = []


    
    for parent, child in tree.edges:
        parent_idx = node_to_idx[parent]
        child_idx = node_to_idx[child]

        length = tree.branch_lengths[(parent, child)]

        #parent to child
        sources.append(parent_idx)
        targets.append(child_idx)
        edge_types.append(parent_to_child)
        edge_attr.append(length)

        sources.append(child_idx)
        targets.append(parent_idx)
        edge_types.append(child_to_parent)
        edge_attr.append(length)

    edge_index = torch.tensor([sources, targets], dtype=torch.long)

    edge_types_tensor = torch.tensor(edge_types, dtype=torch.long)

    edge_attr_tensor = torch.tensor(edge_attr, dtype=torch.float32).unsqueeze(-1)

    return edge_index, edge_types_tensor, edge_attr_tensor









