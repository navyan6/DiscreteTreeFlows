"""
Phase 3: Reference evolutionary process P^0.

Implements mutation prior Q^0 (via ``src/r0_backends``), fitness-biased mutations
Q^0_F (§4.2 Option A site_local or Option B full_esm), and Poisson branching
intensity λ(x) (Alg. 3).

Fitness tilting used by train/gen lives in ``src/bridge/fitness_tilt.py``.
Multi-pLM / substitution R0 wrappers live in ``src/r0_backends.py``.

TreeSBM train/gen still use learned RateHeads for branching by default; set
``branching_mode='poisson_ref'`` on reference rollouts (or gen CLI) to use
Poisson(λ Δt) under P^0 as in paper Alg. 3 / Table D.2.
"""
from __future__ import annotations

from typing import Optional, Protocol, runtime_checkable

import numpy as np
import torch
import torch.nn as nn

from src.tree_state import TreeState
from src.bridge.fitness_tilt import (
    TILT_FULL_ESM,
    TILT_SITE_LOCAL,
    make_sequence_pll_scorer,
    tilt_log_R0_by_fitness,
)
from src.r0_backends import AA_VOCAB, R0Backend


@runtime_checkable
class BranchingIntensity(Protocol):
    def __call__(self, seq: str) -> float: ...


class BranchingIntensityMLP(nn.Module):
    """Single-layer MLP: sequence embedding → branching intensity λ > 0."""

    def __init__(self, embed_dim: int = 1024, esm_c_dim: Optional[int] = None):
        # ``esm_c_dim`` kept as alias for older call sites / tests.
        super().__init__()
        if esm_c_dim is not None:
            embed_dim = int(esm_c_dim)
        self.linear = nn.Linear(embed_dim, 1)

    def forward(self, embeddings: torch.Tensor) -> torch.Tensor:
        """
        Args:
            embeddings: (batch, embed_dim)

        Returns:
            lambdas: (batch,) positive branching rates
        """
        logits = self.linear(embeddings)
        return torch.nn.functional.softplus(logits).squeeze(-1)


class ConstantBranchingIntensity:
    """Sequence-independent Poisson rate λ (D.2 'no seq.-dep. branching' prior)."""

    def __init__(self, lam: float = 1.0):
        self.lam = float(lam)

    def __call__(self, seq: str) -> float:
        return max(self.lam, 1e-6)


class MLPBranchingIntensity:
    """Wrap a trained BranchingIntensityMLP + embedding fn as λ(x)."""

    def __init__(self, mlp: BranchingIntensityMLP, embed_fn):
        self.mlp = mlp
        self.embed_fn = embed_fn
        self.mlp.eval()

    def __call__(self, seq: str) -> float:
        with torch.no_grad():
            emb = self.embed_fn(seq)
            if emb.dim() == 1:
                emb = emb.unsqueeze(0)
            return max(float(self.mlp(emb)[0].item()), 1e-6)


def sample_poisson_offspring(lam: float, dt: float, rng: Optional[np.random.Generator] = None) -> int:
    """b ~ Poisson(λ Δt) as in Alg. 3."""
    mean = max(float(lam), 0.0) * float(dt)
    if rng is None:
        return int(np.random.poisson(mean))
    return int(rng.poisson(mean))


class ReferenceProcess:
    """
    Biologically grounded reference process P^0(T_{0:1} | x0).

    - Mutation prior Q^0 from any ``R0Backend`` (ESM-2 / ESM-C / JTT / …)
    - Fitness tilt Q^0_F via Option A (``site_local``) or Option B (``full_esm``)
    - Branching b ~ Poisson(λ(x) Δt) via ``branching_intensity``
    """

    def __init__(
        self,
        r0_backend: R0Backend,
        branching_intensity: Optional[BranchingIntensity] = None,
        beta: float = 1.0,
        fitness_score: str = "log_R0",
        fitness_tilt_mode: str = TILT_SITE_LOCAL,
        fitness_esm_batch_size: int = 8,
        fitness_esm_top_k: Optional[int] = None,
        p_stop: float = 0.0,
        rng: Optional[np.random.Generator] = None,
        # Legacy kwargs kept so older call sites / docs still construct:
        esm2_model=None,
        esm2_alphabet=None,
        esm_c_model=None,
        branching_mlp: Optional[BranchingIntensityMLP] = None,
    ):
        if r0_backend is None and esm2_model is not None:
            raise ValueError(
                "Legacy esm2_model/esm2_alphabet ctor is removed; "
                "pass an R0Backend from src.r0_backends.build_r0_backend(...)."
            )
        self.r0 = r0_backend
        self.beta = float(beta)
        self.fitness_score = fitness_score
        self.fitness_tilt_mode = fitness_tilt_mode
        self.fitness_esm_batch_size = int(fitness_esm_batch_size)
        self.fitness_esm_top_k = fitness_esm_top_k
        self.p_stop = float(p_stop)
        self.rng = rng if rng is not None else np.random.default_rng()
        self._fitness_cache: dict = {}
        self._fitness_scorer = None
        if self.fitness_tilt_mode == TILT_FULL_ESM and self.beta != 0.0:
            self._fitness_scorer = make_sequence_pll_scorer(self.r0)

        if branching_intensity is not None:
            self.branching_intensity = branching_intensity
        elif branching_mlp is not None and esm_c_model is not None:
            # Legacy path: MLP on ESM-C token embedding.
            def _embed(seq: str) -> torch.Tensor:
                tokens = esm_c_model.tokenize_seq(seq)
                tokens = torch.tensor(tokens).unsqueeze(0)
                return esm_c_model.encode(tokens)[0, 0, :]

            self.branching_intensity = MLPBranchingIntensity(branching_mlp, _embed)
        else:
            self.branching_intensity = ConstantBranchingIntensity(1.0)

    def get_mutation_log_rates(self, seq: str) -> torch.Tensor:
        """Return tilted log R0 for one sequence: [L, 20]."""
        log_R0 = self.r0.log_mutation_rates([seq], max_seq_len=len(seq))  # [1, L, 20]
        return tilt_log_R0_by_fitness(
            log_R0,
            beta=self.beta,
            score=self.fitness_score,
            mode=self.fitness_tilt_mode,
            sequences=[seq] if self.fitness_tilt_mode == TILT_FULL_ESM else None,
            fitness_scorer=self._fitness_scorer,
            cache=self._fitness_cache,
            batch_size=self.fitness_esm_batch_size,
            top_k_aas=self.fitness_esm_top_k,
        )[0]

    def get_mutation_rates(self, seq: str) -> np.ndarray:
        """Q^0_F probabilities [L, 20]."""
        return self.get_mutation_log_rates(seq).exp().cpu().numpy()

    def get_branching_intensity(self, seq: str) -> float:
        return max(float(self.branching_intensity(seq)), 1e-6)

    def rollout(
        self,
        x0: str,
        horizon: float = 1.0,
        dt: float = 0.05,
        max_nodes: int = 256,
    ) -> TreeState:
        """
        Algorithm 3: ReferenceRollout.
        Mutate via tilted pLM/substitution prior; branch via Poisson(λ Δt).
        """
        tree = TreeState.root_only(x0)

        for _t in np.arange(0, horizon, dt):
            if len(tree.node_ids) >= max_nodes:
                break

            new_active_leaves: list[str] = []

            for leaf_id in list(tree.active_leaves):
                leaf_seq = tree.node_seqs[leaf_id]

                # Optional lineage stop (Alg. 3 line 12).
                if self.p_stop > 0.0 and self.rng.random() < self.p_stop * dt:
                    continue

                lam = self.get_branching_intensity(leaf_seq)
                b = sample_poisson_offspring(lam, dt, self.rng)

                if b > 0:
                    child_seqs = [self._mutate_sequence(leaf_seq) for _ in range(b)]
                    tree = tree.branch_node(leaf_id, child_seqs)
                    new_active_leaves.extend(
                        [f"{leaf_id}_child_{i}" for i in range(b)]
                    )
                else:
                    # No branching: mutate in place (local edit along lineage).
                    mut_seq = self._mutate_sequence(leaf_seq)
                    tree = self._replace_seq(tree, leaf_id, mut_seq)
                    new_active_leaves.append(leaf_id)

            for leaf_id in new_active_leaves:
                if leaf_id == tree.root_id or tree.get_parent(leaf_id) is None:
                    continue
                tree = tree.extend_branch(leaf_id, dt)

            tree = TreeState(
                node_ids=tree.node_ids,
                root_id=tree.root_id,
                edges=tree.edges,
                branch_lengths=tree.branch_lengths,
                node_seqs=tree.node_seqs,
                active_leaves=new_active_leaves,
            )

        return tree

    @staticmethod
    def _replace_seq(tree: TreeState, node_id: str, new_seq: str) -> TreeState:
        seqs = dict(tree.node_seqs)
        seqs[node_id] = new_seq
        return TreeState(
            node_ids=tree.node_ids.copy(),
            root_id=tree.root_id,
            edges=tree.edges.copy(),
            branch_lengths=tree.branch_lengths.copy(),
            node_seqs=seqs,
            active_leaves=tree.active_leaves.copy(),
        )

    def _mutate_sequence(self, seq: str) -> str:
        """Sample one edit from tilted Q^0_F site distribution."""
        qf_probs = self.get_mutation_rates(seq)  # (L, 20)
        L = min(len(seq), qf_probs.shape[0])
        if L == 0:
            return seq

        pos_probs = qf_probs[:L].sum(axis=1)
        pos_probs = np.clip(pos_probs, 0, None)
        if pos_probs.sum() <= 0:
            return seq
        pos_probs /= pos_probs.sum()
        pos = int(self.rng.choice(L, p=pos_probs))

        target_probs = np.clip(qf_probs[pos], 0, None)
        if target_probs.sum() <= 0:
            return seq
        target_probs /= target_probs.sum()
        target_aa = str(self.rng.choice(list(AA_VOCAB), p=target_probs))
        return seq[:pos] + target_aa + seq[pos + 1 :]

    def _get_fitness(self, seq: str) -> float:
        """Mean site log-prob under untilted R0 (Option B / diagnostics)."""
        log_R0 = self.r0.log_mutation_rates([seq], max_seq_len=len(seq))[0]
        return float(log_R0.mean().item())
