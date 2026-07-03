#causal masking for tree transformer

import torch

def build_causal_attention_mask(
    node_ids: list[str],
    node_times: dict[str, float],
) -> torch.Tensor:
    """
    Temporal causal mask: node i attends to node j only if time(j) <= time(i).
    """
    times = torch.tensor([node_times[nid] for nid in node_ids], dtype=torch.float32)
    mask = times.unsqueeze(1) >= times.unsqueeze(0)  # [N, N]
    return mask

# Alias used in networks.py
build_temporal_attention_mask = build_causal_attention_mask

#based on the treetime node id

