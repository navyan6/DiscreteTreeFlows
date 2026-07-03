from dataclasses import dataclass
import torch


@dataclass
class TreeEncoderInput:
    """Tensor representation of a tree state T = (V, E, r, l, X)."""

    node_ids: list[str]
    x: torch.Tensor
    structural_features: torch.Tensor
    lap_pe: torch.Tensor
    edge_index: torch.Tensor
    edge_type: torch.Tensor
    edge_attr: torch.Tensor
    root_index: int


    