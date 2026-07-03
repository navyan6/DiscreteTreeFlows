import torch

from src.tree_state import TreeState


#compute the max depth of a leave from a root using a stack (dfs)
def compute_depth(tree: TreeState) -> dict[str, int]:
    depths = {tree.root_id: 0}
    stack = [tree.root_id]
    while stack:
        current_node = stack.pop()

        for child in tree.get_children(current_node):
            depths[child] = depths[current_node] + 1
            stack.append(child)

    if len(depths) != tree.n_nodes():
        raise ValueError("Some nodes are not reachable by root")


    return depths


def compute_structural_features(
        tree: TreeState,
        node_to_index: dict[str, int],
) -> torch.Tensor:

    """return: normalized depth, root indicator, leaf indicator
    """

    depths = compute_depth(tree)
    max_depth = max(depths.values(), default=0)

    #create a matrix of storing one structural feature per node 

    features = torch.zeros(tree.n_nodes(), 3, dtype=torch.float32)

    for node_id in tree.node_ids: 
        index = node_to_index[node_id]

        normalized_depth = (
            depths[node_id] / max_depth
            if max_depth > 0
            else 0.0
        )

        features[index, 0] = normalized_depth
        features[index, 1] = float(node_id == tree.root_id)
        features[index, 2] = float(tree.is_leaf(node_id))

    return features





        


