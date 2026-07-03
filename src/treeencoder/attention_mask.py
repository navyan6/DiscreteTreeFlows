#causal masking for tree transformer

import torch


def causal_mask(seq_len, device):
    mask = torch.tril(torch.ones((seq_len, seq_len), device=device)).unsqueeze(0)
    return mask


def build_temporal_attention_mask(
    node_ids: list[str],
    node_times: dict[str, float],
) -> torch.Tensor:
    """
    Temporal causal mask: node i attends to node j only if time(j) <= time(i).

    Args:
        node_ids:   ordered list of node IDs (defines row/col ordering)
        node_times: dict mapping node_id -> numdate (decimal year)

    Returns:
        [N, N] bool tensor — mask[i, j] = True means i can attend to j
    """
    times = torch.tensor([node_times[nid] for nid in node_ids], dtype=torch.float32)
    # mask[i, j] = (times[j] <= times[i])
    mask = times.unsqueeze(1) >= times.unsqueeze(0)  # [N, N]
    return mask



