import torch

from src.tree_state import TreeState

def compute_laplacian_pe(
    tree: TreeState,
    node_to_index: dict[str, int],
    num_eigenvectors: int,
    device: str | torch.device = "cpu",
) -> torch.Tensor:
    """
    Compute Laplacian positional encodings for every node.
    """
    if num_eigenvectors < 0:
        raise ValueError("num_eigenvectors must be nonnegative")

    num_nodes = tree.n_nodes()

    if len(node_to_index) != num_nodes:
        raise ValueError(
            "node_to_index must contain one entry for every tree node"
        )

    if set(node_to_index) != set(tree.node_ids):
        raise ValueError(
            "node_to_index keys must exactly match tree.node_ids"
        )
    if num_eigenvectors == 0:
        return torch.empty(
            (num_nodes, 0),
            dtype=torch.float32,
        )
    if num_nodes == 1:
        return torch.zeros(
            (1, num_eigenvectors),
            dtype=torch.float32,
        )

    # A[i, j] = 1 when nodes i and j share a tree edge.
    adjacency = torch.zeros(
        (num_nodes, num_nodes),
        dtype=torch.float32,
    )

    for parent_id, child_id in tree.edges:
        parent_index = node_to_index[parent_id]
        child_index = node_to_index[child_id]

        # Treat the directed tree edge as undirected for the Laplacian.
        adjacency[parent_index, child_index] = 1.0
        adjacency[child_index, parent_index] = 1.0

    # Degree of node i = number of neighbors of node i.
    degree = adjacency.sum(dim=1)

    inverse_sqrt_degree = torch.zeros_like(degree)
    nonzero_degree = degree > 0
    inverse_sqrt_degree[nonzero_degree] = (
        degree[nonzero_degree].pow(-0.5)
    )

    degree_inv_sqrt = torch.diag(inverse_sqrt_degree)
    identity = torch.eye(num_nodes, dtype=torch.float32)

    laplacian = (
        identity
        - degree_inv_sqrt
        @ adjacency
        @ degree_inv_sqrt
    )

    # Eigendecomp on CPU — GPU eigh is numerically unstable on B200 (sm_100)
    eigenvalues, eigenvectors = torch.linalg.eigh(laplacian)

    available = min(
        num_eigenvectors,
        num_nodes - 1,
    )

    positional_encoding = eigenvectors[
        :,
        1 : 1 + available,
    ]

    if available < num_eigenvectors:
        padding = torch.zeros(
            (num_nodes, num_eigenvectors - available),
            dtype=positional_encoding.dtype,
        )

        positional_encoding = torch.cat(
            [positional_encoding, padding],
            dim=1,
        )

    return positional_encoding.to(device)