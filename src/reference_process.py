"""
Phase 3: Reference evolutionary process P^0.
Implements mutation prior Q^0, fitness-biased mutations Q^0_F, and branching intensity λ(x).
"""
import torch
import torch.nn as nn
import numpy as np
from typing import Optional

from src.tree_state import TreeState


class BranchingIntensityMLP(nn.Module):
    """Single-layer MLP: ESM-C embedding gives branching intensity λ."""

    def __init__(self, esm_c_dim: int = 1024):
        super().__init__()
        self.linear = nn.Linear(esm_c_dim, 1)

    def forward(self, embeddings: torch.Tensor) -> torch.Tensor:
        """
        Args:
            embeddings: (batch, esm_c_dim) ESM-C embeddings

        Returns:
            lambdas: (batch,) positive branching rates
        """
        logits = self.linear(embeddings)  # (batch, 1)
        lambdas = torch.nn.functional.softplus(logits).squeeze(-1)  # (batch,)
        return lambdas


class ReferenceProcess:
    """
    Biologically grounded reference process P^0(T_{0:1} | x0).
    Combines:
    - Mutation prior Q^0(x, x') from ESM-2 masked LM
    - Fitness-biased mutations Q^0_F(x, x') ∝ Q^0(x, x') exp(βF(x'))
    - Branching intensity λ(x_v) from trained MLP
    """

    def __init__(
        self,
        esm2_model,
        esm2_alphabet,
        esm_c_model,
        branching_mlp: BranchingIntensityMLP,
        beta: float = 1.0,
    ):
        """
        Args:
            esm2_model: Loaded ESM-2 model for mutations
            esm2_alphabet: ESM-2 alphabet for tokenization
            esm_c_model: Loaded ESM-C model for embeddings
            branching_mlp: Trained 1-layer MLP for λ(x)
            beta: Fitness temperature parameter
        """
        self.esm2_model = esm2_model
        self.esm2_alphabet = esm2_alphabet
        self.esm_c_model = esm_c_model
        self.branching_mlp = branching_mlp
        self.beta = beta

        # Freeze all models
        for model in [esm2_model, esm_c_model]:
            model.eval()
            for param in model.parameters():
                param.requires_grad = False

    def get_mutation_rates(self, seq: str) -> np.ndarray:
        """
        Compute Q^0_F(x, x') for all single-residue mutations.

        Returns:
            mutation_rates: (L, 20) array of transition rates
                where mutation_rates[i, j] = rate of position i → amino acid j
        """
        with torch.no_grad():
            # Tokenize
            tokens = self.esm2_alphabet.encode(seq)
            tokens = torch.tensor(tokens).unsqueeze(0)  # (1, L)

            # Get logits from ESM-2
            results = self.esm2_model(tokens, repr_layers=[33])
            logits = results["logits"][0]  # (L, 33 vocab)

            # Extract amino acid logits (indices 4-23 for 20 standard AAs)
            aa_logits = logits[:, 4:24]  # (L, 20)

            # Convert to probabilities
            aa_probs = torch.softmax(aa_logits, dim=-1)  # (L, 20)

            # Get fitness for this sequence: mean per-position log prob
            aa_logprobs = torch.log_softmax(aa_logits, dim=-1)
            fitness = aa_logprobs.mean().item()  # scalar

            # Bias mutations by fitness of target: Q^0_F(x, x') ∝ Q^0(x, x') exp(β F(x'))
            # For simplicity: use current sequence fitness as proxy
            # Better: would compute F(x') for each mutation, but expensive
            fitness_factor = np.exp(self.beta * fitness)

            mutation_rates = (aa_probs.cpu().numpy() * fitness_factor)

        return mutation_rates  # (L, 20)

    def get_branching_intensity(self, seq: str) -> float:
        """
        Get branching intensity λ(x_v) for a sequence.

        Args:
            seq: amino acid sequence

        Returns:
            lambda: positive scalar, expected offspring per unit time
        """
        with torch.no_grad():
            # Embed with ESM-C
            tokens = self.esm_c_model.tokenize_seq(seq)
            tokens = torch.tensor(tokens).unsqueeze(0)  # (1, L)
            embedding = self.esm_c_model.encode(tokens)[0, 0, :]  # (1024,) - use first token

            # Pass through MLP
            embedding = embedding.unsqueeze(0)  # (1, 1024)
            lambda_val = self.branching_mlp(embedding).item()

        return max(lambda_val, 1e-6)  # ensure positive

    def rollout(
        self, x0: str, horizon: float = 1.0, dt: float = 0.05, max_nodes: int = 256
    ) -> TreeState:
        """
        Algorithm 3: ReferenceRollout.
        Generate tree from root sequence under P^0.

        Args:
            x0: root sequence
            horizon: maximum time
            dt: time step size
            max_nodes: stop if tree grows too large

        Returns:
            tree: final TreeState after rollout
        """
        tree = TreeState.root_only(x0)
        node_counter = 0

        for t in np.arange(0, horizon, dt):
            if len(tree.node_ids) >= max_nodes:
                break

            new_active_leaves = []

            for leaf_id in tree.active_leaves:
                leaf_seq = tree.node_seqs[leaf_id]

                # Get branching rate for this sequence
                lambda_leaf = self.get_branching_intensity(leaf_seq)

                # Sample offspring count from Poisson(λ * dt)
                b = np.random.poisson(lambda_leaf * dt)

                if b > 0:
                    # Generate offspring with mutations (Algorithm 3, Lines 7-8)
                    child_seqs = []
                    for _ in range(b):
                        child_seq = self._mutate_sequence(leaf_seq)
                        child_seqs.append(child_seq)

                    # Add to tree
                    tree = tree.branch_node(leaf_id, child_seqs)
                    new_active_leaves.extend(
                        [f"{leaf_id}_child_{i}" for i in range(b)]
                    )
                else:
                    # No branching, keep growing this leaf
                    new_active_leaves.append(leaf_id)

            # Algorithm 3, Line 11: Update existing branch lengths by Δt
            for leaf_id in new_active_leaves:
                tree = tree.extend_branch(leaf_id, dt)

            tree.active_leaves = new_active_leaves

        return tree

    def _mutate_sequence(self, seq: str) -> str:
        """
        Algorithm 3, Lines 7-8: Sample mutation and accept/reject by fitness.

        1. Sample edit proposal x' from p_pLM (Q^0)
        2. Accept x' with probability ∝ exp(βF(x'))

        Args:
            seq: current sequence

        Returns:
            mutated sequence (after accept/reject)
        """
        aa_alphabet = "ACDEFGHIKLMNPQRSTVWY"

        with torch.no_grad():
            # Get Q^0 (ESM-2 mutation probabilities)
            tokens = self.esm2_alphabet.encode(seq)
            tokens = torch.tensor(tokens).unsqueeze(0)
            results = self.esm2_model(tokens, repr_layers=[33])
            logits = results["logits"][0]
            aa_logits = logits[:, 4:24]
            q0_probs = torch.softmax(aa_logits, dim=-1).cpu().numpy()  # (L, 20)

        # Sample position and target amino acid from Q^0
        pos_probs = q0_probs.sum(axis=1)
        pos_probs /= pos_probs.sum()
        pos = np.random.choice(len(seq), p=pos_probs)

        target_probs = q0_probs[pos]
        target_probs /= target_probs.sum()
        target_aa = np.random.choice(aa_alphabet, p=target_probs)

        # Create mutant
        mutant = seq[:pos] + target_aa + seq[pos+1:]

        # Accept/reject: compute F(x') and accept with prob ∝ exp(βF(x'))
        fitness_mutant = self._get_fitness(mutant)
        accept_prob = np.exp(self.beta * fitness_mutant)

        # Cap acceptance probability at 1
        accept_prob = min(accept_prob, 1.0)

        if np.random.rand() < accept_prob:
            return mutant
        else:
            return seq  # Reject: stay at current sequence

    def _get_fitness(self, seq: str) -> float:
        """
        Compute fitness F(x) = mean per-position log-probability under ESM-2.

        Args:
            seq: amino acid sequence

        Returns:
            fitness: scalar (can be negative)
        """
        with torch.no_grad():
            tokens = self.esm2_alphabet.encode(seq)
            tokens = torch.tensor(tokens).unsqueeze(0)
            results = self.esm2_model(tokens, repr_layers=[33])
            logits = results["logits"][0]
            aa_logits = logits[:, 4:24]
            aa_logprobs = torch.log_softmax(aa_logits, dim=-1)
            fitness = aa_logprobs.mean().item()
        return fitness
