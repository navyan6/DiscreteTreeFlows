"""
Encodes tree state and predicts controlled rates for bridge matching.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Tuple, Optional
from src.treeencoder.attention_mask import build_temporal_attention_mask


class EdgeWeightingMLP(nn.Module):
    """Learn per-head edge weights from branch lengths."""

    def __init__(self, d_model: int, n_heads: int):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(1, 64),
            nn.ReLU(),
            nn.Linear(64, n_heads),
            nn.Sigmoid(),
        )

    def forward(self, branch_lengths: torch.Tensor) -> torch.Tensor:

        return self.mlp(branch_lengths.unsqueeze(-1))


def scaled_dot_product_attention_with_edges(
    q: torch.Tensor,
    k: torch.Tensor,
    v: torch.Tensor,
    attn_mask: torch.Tensor,
    edge_index: Optional[torch.Tensor] = None,
    edge_weights_per_head: Optional[torch.Tensor] = None,
    dropout_p: float = 0.0,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Scaled dot-product attention with sparse mask and learned edge weighting.
    """
    d_k = q.shape[-1]
    batch_size, n_heads, seq_len, _ = q.shape
    scores = torch.matmul(q, k.transpose(-2, -1)) / (d_k ** 0.5)

    # Apply sparse mask: set masked positions to -inf
    scores = scores.masked_fill(~attn_mask.unsqueeze(0).unsqueeze(0), -1e9)

    # Apply learned edge weights if provided
    if edge_index is not None and edge_weights_per_head is not None:
        # Build weight matrix [seq_len, seq_len, n_heads]
        edge_weight_matrix = torch.ones(
            (seq_len, seq_len, n_heads),
            device=scores.device,
            dtype=scores.dtype,
        )

        # Apply weights to edge positions
        src_idx, tgt_idx = edge_index[0], edge_index[1]
        for e in range(src_idx.shape[0]):
            i, j = src_idx[e].item(), tgt_idx[e].item()
            edge_weight_matrix[i, j, :] = edge_weights_per_head[e, :]

        # Modulate scores: [batch, n_heads, seq_len, seq_len] *= [seq_len, seq_len, n_heads]
        # Need to transpose to match dimensions
        edge_weight_matrix = edge_weight_matrix.permute(2, 0, 1)  # [n_heads, seq_len, seq_len]
        scores = scores * edge_weight_matrix.unsqueeze(0)  # broadcast batch dimension

    attn_weights = torch.softmax(scores, dim=-1)
    attn_weights = torch.dropout(attn_weights, dropout_p, train=True)

    attn_output = torch.matmul(attn_weights, v)
    return attn_output, attn_weights


class TreeEncoder(nn.Module):
    """
    Graph transformer over a phylogenetic tree.

    Expects pre-computed node embeddings from NodeEncoder
    (ESM2 + structural features + Laplacian PE already fused → d_model).

    Each TransformerGraphLayer runs multi-head attention masked by a temporal
    causal mask (node i attends to node j only if time(j) <= time(i)), with
    edge weights learned from branch lengths via EdgeWeightingMLP. Followed by
    a position-wise MLP with residual connections — standard transformer block.
    """

    def __init__(
        self,
        d_model: int = 128,
        n_layers: int = 4,
        n_heads: int = 8,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.d_model = d_model

        self.transformer_layers = nn.ModuleList(
            [TransformerGraphLayer(d_model, n_heads, dropout) for _ in range(n_layers)]
        )

        self.graph_readout = nn.Sequential(
            nn.Linear(d_model, d_model),
            nn.ReLU(),
            nn.Linear(d_model, d_model),
        )

        # Time conditioning: inject interpolation time t into node features
        self.time_proj = nn.Sequential(
            nn.Linear(1, d_model),
            nn.SiLU(),
            nn.Linear(d_model, d_model),
        )

    def forward(
        self,
        h: torch.Tensor,
        node_ids: list[str],
        node_times: dict[str, float],
        edge_index: torch.Tensor,
        branch_lengths: torch.Tensor,
        t_scalar: Optional[float] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Args:
            h:              [N, d_model] pre-computed node embeddings (from NodeEncoder)
            node_ids:       ordered list of node IDs (length N)
            node_times:     dict node_id -> numdate (decimal year) for temporal causal mask
            edge_index:     [2, E] edge indices (both directions from build_edges)
            branch_lengths: [E] branch length per edge for EdgeWeightingMLP
            t_scalar:       interpolation time in [0,1]; None skips time conditioning

        Returns:
            H_T:     [N, d_model] contextualized node embeddings
            h_graph: [d_model]   mean-pooled graph embedding
        """
        attn_mask = build_temporal_attention_mask(node_ids, node_times).to(h.device)

        # Inject bridge time into every node embedding before attention
        if t_scalar is not None:
            t_in = torch.tensor([[t_scalar]], dtype=h.dtype, device=h.device)
            t_emb = self.time_proj(t_in)  # [1, d_model]
            h = h + t_emb                 # broadcast over N nodes

        for layer in self.transformer_layers:
            h = layer(h, attn_mask, edge_index, branch_lengths)

        h_graph = self.graph_readout(h.mean(dim=0))

        return h, h_graph


class TransformerGraphLayer(nn.Module):
    """Single transformer layer with sparse graph attention and learned edge weighting."""

    def __init__(self, d_model: int = 256, n_heads: int = 8, dropout: float = 0.1):
        super().__init__()
        self.d_model = d_model
        self.n_heads = n_heads
        self.dropout_p = dropout

        self.q_proj = nn.Linear(d_model, d_model)
        self.k_proj = nn.Linear(d_model, d_model)
        self.v_proj = nn.Linear(d_model, d_model)
        self.out_proj = nn.Linear(d_model, d_model)

        self.edge_weighting = EdgeWeightingMLP(d_model, n_heads)

        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)

        self.mlp = nn.Sequential(
            nn.Linear(d_model, 4 * d_model),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(4 * d_model, d_model),
        )

    def forward(
        self,
        h: torch.Tensor,
        attn_mask: torch.Tensor,
        edge_index: Optional[torch.Tensor] = None,
        branch_lengths: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:

        # Project to Q, K, V
        q = self.q_proj(h)  # (n_nodes, d_model)
        k = self.k_proj(h)
        v = self.v_proj(h)

        # Reshape for multi-head attention
        batch_size = 1
        seq_len = h.shape[0]
        q = q.view(batch_size, seq_len, self.n_heads, -1).transpose(1, 2)
        k = k.view(batch_size, seq_len, self.n_heads, -1).transpose(1, 2)
        v = v.view(batch_size, seq_len, self.n_heads, -1).transpose(1, 2)

        # Compute learned edge weights
        edge_weights_per_head = None
        if edge_index is not None and branch_lengths is not None:
            edge_weights_per_head = self.edge_weighting(branch_lengths)  # (n_edges, n_heads)

        # Sparse attention with mask and learned edge weights
        attn_out, _ = scaled_dot_product_attention_with_edges(
            q,
            k,
            v,
            attn_mask,
            edge_index=edge_index,
            edge_weights_per_head=edge_weights_per_head,
            dropout_p=self.dropout_p,
        )

        # Reshape back
        attn_out = attn_out.transpose(1, 2).contiguous()
        attn_out = attn_out.view(batch_size, seq_len, self.d_model).squeeze(0)

        attn_out = self.out_proj(attn_out)

        # Residual + layer norm
        h = self.norm1(h + attn_out)

        # MLP
        mlp_output = self.mlp(h)
        h = self.norm2(h + mlp_output)

        return h


class RateHeads(nn.Module):
    """
    Four prediction heads for controlled rates R_θ(T, T', t).

    outputs:
    - Mutation rates: (n_active, L, 20) logits over position × amino acid
    - Branching rates: (n_active,) Poisson mean λ
    - Branch length: (n_active,) continuous extension magnitude
    - Stop probability: (n_active,) Bernoulli probability
    """

    def __init__(self, d_model: int = 256, max_seq_len: int = 512):
        super().__init__()
        self.d_model = d_model
        self.max_seq_len = max_seq_len
        self.mutation_head = nn.Sequential(
            nn.Linear(d_model, d_model),
            nn.ReLU(),
            nn.Linear(d_model, max_seq_len * 20),
        )

        # Branching head: Poisson parameter λ
        self.branching_head = nn.Sequential(
            nn.Linear(d_model, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Softplus(),  # ensures λ > 0
        )

        # Branch length head: continuous magnitude
        self.branch_length_head = nn.Sequential(
            nn.Linear(d_model, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Softplus(),
        )

        # Stop head: stopping probability
        self.stop_head = nn.Sequential(
            nn.Linear(d_model, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid(),
        )

    def forward(
        self, H_T: torch.Tensor, active_leaf_indices: list[int]
    ) -> dict[str, torch.Tensor]:
        """
        Predict rates for active leaves only.
        """
        # Extract active leaf embeddings
        h_active = H_T[active_leaf_indices]  # (n_active, d_model)

        # Mutation head
        mutation_logits = self.mutation_head(h_active)  # (n_active, L*20)
        mutation_logits = mutation_logits.reshape(
            -1, self.max_seq_len, 20
        )  # (n_active, L, 20)

        # Branching head
        branching_rate = self.branching_head(h_active).squeeze(-1)  # (n_active,)

        # Branch length head
        branch_length = self.branch_length_head(h_active).squeeze(-1)  # (n_active,)

        # Stop head
        stop_prob = self.stop_head(h_active).squeeze(-1)  # (n_active,)

        return {
            "mutation_logits": mutation_logits,
            "branching_rate": branching_rate,
            "branch_length": branch_length,
            "stop_prob": stop_prob,
        }
