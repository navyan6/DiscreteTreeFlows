"""
Node encoder: fuse PLM embeddings, structural features, and Laplacian PE.
"""
import torch
import torch.nn as nn


class NodeEncoder(nn.Module):
    def __init__(
        self,
        d_plm: int = 320,
        d_struct: int = 3,
        d_laplacian: int = 8,
        d_node: int = 128,
        activation: str = "relu",
    ):
        #3 parts + activation function
        super().__init__()
        self.d_plm = d_plm
        self.d_struct = d_struct
        self.d_laplacian = d_laplacian
        self.d_node = d_node

        input_dim = d_plm + d_struct + d_laplacian

        self.projection = nn.Linear(input_dim, d_node)
#relu
        if activation.lower() == "relu":
            self.activation = nn.ReLU()
        elif activation.lower() == "gelu":
            self.activation = nn.GELU()
        elif activation.lower() == "none":
            self.activation = nn.Identity()
        else:
            raise ValueError(f"Unknown activation: {activation}")

    def forward(
        self,
        x: torch.Tensor,
        structural_features: torch.Tensor,
        lap_pe: torch.Tensor,
    ) -> torch.Tensor:
        #embeddings from esmc
        # Validate input shapes
        batch_size = x.shape[0]
        assert structural_features.shape == (batch_size, self.d_struct), (
            f"Expected structural_features shape ({batch_size}, {self.d_struct}), "
            f"got {structural_features.shape}"
        )
        assert lap_pe.shape == (batch_size, self.d_laplacian), (
            f"Expected lap_pe shape ({batch_size}, {self.d_laplacian}), "
            f"got {lap_pe.shape}"
        )

        combined = torch.cat([x, structural_features, lap_pe], dim=-1)

        h = self.projection(combined)

        h = self.activation(h)

        return h
