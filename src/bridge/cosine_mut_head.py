"""
CoSiNE-style site-local mutation rates for TreeSBM RateHeads.

Unlike the residual head (log R_θ = log R0 + c_θ), this CNN predicts
per-site AA logits from parent identity + local k-mer context + IMGT
position + a projected tree vector. A learned scalar gate can mix in
frozen R0; default init makes the mix ≈ 0 so ESM conservation does not
dominate.

Antibody Recipe C only. Viral RateHeads stay residual.
"""

from __future__ import annotations

import math

import torch
import torch.nn as nn


class CosineSiteMutHead(nn.Module):
    def __init__(
        self,
        d_model: int = 128,
        d_aa: int = 32,
        d_pos: int = 32,
        n_channels: int = 64,
        kernel_size: int = 5,
        n_layers: int = 2,
        max_seq_len: int = 160,
        r0_mix_init: float = -4.0,
    ):
        super().__init__()
        if kernel_size % 2 != 1:
            raise ValueError("kernel_size must be odd")
        self.max_seq_len = int(max_seq_len)
        self.aa_emb = nn.Embedding(21, d_aa)
        self.pos_emb = nn.Embedding(self.max_seq_len, d_pos)
        self.tree_proj = nn.Linear(d_model, d_aa)
        in_ch = d_aa + d_pos + d_aa
        layers: list[nn.Module] = []
        ch = in_ch
        pad = kernel_size // 2
        for _ in range(n_layers):
            layers.append(nn.Conv1d(ch, n_channels, kernel_size, padding=pad))
            layers.append(nn.GELU())
            ch = n_channels
        self.conv = nn.Sequential(*layers)
        self.out = nn.Conv1d(n_channels, 20, kernel_size=1)
        # sigmoid(r0_mix_logit) weight on frozen R0. -4 → ~0.018.
        self.r0_mix_logit = nn.Parameter(torch.tensor(float(r0_mix_init)))

    def forward(
        self,
        h_active: torch.Tensor,
        aa_indices: torch.Tensor,
        log_R0_mut: torch.Tensor,
    ) -> torch.Tensor:
        n, L, _ = log_R0_mut.shape
        if aa_indices is None:
            raise ValueError("CosineSiteMutHead requires aa_indices [n, L]")
        aa = aa_indices.to(device=h_active.device, dtype=torch.long)
        if aa.shape != (n, L):
            raise ValueError(f"aa_indices shape {tuple(aa.shape)} != {(n, L)}")
        if L > self.max_seq_len:
            raise ValueError(f"L={L} exceeds max_seq_len={self.max_seq_len}")
        aa_h = self.aa_emb(aa.clamp(0, 20))
        pos_ids = torch.arange(L, device=h_active.device)
        pos_h = self.pos_emb(pos_ids).unsqueeze(0).expand(n, -1, -1)
        tree_h = self.tree_proj(h_active).unsqueeze(1).expand(-1, L, -1)
        x = torch.cat([aa_h, pos_h, tree_h], dim=-1)  # [n, L, C]
        x = x.transpose(1, 2)  # [n, C, L]
        logits = self.out(self.conv(x)).transpose(1, 2)  # [n, L, 20]
        mix = torch.sigmoid(self.r0_mix_logit)
        return (1.0 - mix) * logits + mix * log_R0_mut


def cosine_r0_mix(module: CosineSiteMutHead) -> float:
    return float(torch.sigmoid(module.r0_mix_logit).detach().cpu())


def default_r0_mix_init() -> float:
    return -4.0


def expected_init_mix(logit: float = -4.0) -> float:
    return 1.0 / (1.0 + math.exp(-logit))
