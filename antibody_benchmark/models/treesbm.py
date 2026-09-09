"""TreeSBM adapter: sequence evolution on FIXED observed topology+BL."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

from antibody_benchmark.models.base import EvolutionModel


class TreeSBMModel(EvolutionModel):
    """Sequence-only comparison: topology and branch lengths are supplied externally.

    ``sample_child`` evolves parent→child for the observed branch length using the
    project's RateHeads + mutate_sequence_independent path (no tree generation).
    """

    name = "treesbm"
    alphabet = "aa"

    def __init__(
        self,
        checkpoint: Optional[str] = None,
        force_observed_topology: bool = True,
        device: str | None = None,
        max_seq_len: int = 566,
        mutation_rate_scale: float = 1.0,
    ):
        if not force_observed_topology:
            raise RuntimeError(
                "FAIRNESS VIOLATION: TreeSBM must use force_observed_topology=True "
                "for the shared sequence rollout benchmark."
            )
        self.checkpoint = checkpoint
        self.device = device or "cpu"
        self.max_seq_len = int(max_seq_len)
        self.mutation_rate_scale = float(mutation_rate_scale)
        self._ready = False
        self._load()

    def _load(self) -> None:
        if not self.checkpoint:
            raise RuntimeError(
                "BLOCKER: TreeSBM checkpoint path not set in configs "
                "(models.treesbm.checkpoint)."
            )
        ckpt = Path(self.checkpoint)
        if not ckpt.is_absolute():
            repo_root = Path(__file__).resolve().parents[2]
            ckpt = repo_root / ckpt
        if not ckpt.is_file():
            raise RuntimeError(f"BLOCKER: TreeSBM checkpoint not found: {ckpt}")
        self.checkpoint = str(ckpt)

        repo_root = Path(__file__).resolve().parents[2]
        if str(repo_root) not in sys.path:
            sys.path.insert(0, str(repo_root))

        try:
            import torch
            from transformers import AutoTokenizer, EsmForMaskedLM

            from scripts.eval_single_tree import (
                AA_VOCAB,
                get_lm_logits,
                load_models,
            )
            from src.bridge.fitness_tilt import tilt_log_R0_by_fitness
            from src.bridge.shm_site_prior import apply_shm_site_prior
            from src.bridge.losses import _build_seq_indices
            from src.bridge.mutation_sample import mutate_sequence_independent
            from src.networks import RateHeads  # noqa: F401
            from src.r0_backends import (
                BACKEND_THRIFTY_AA,
                build_r0_backend,
                normalize_backend_name,
            )
            from src.tree_state import TreeState
            from src.treeencoder.edges import build_edges
            from src.treeencoder.laplacian import compute_laplacian_pe
            from src.treeencoder.plm_embeddings import ESM2Embedder
            from src.treeencoder.structural_features import compute_structural_features
        except Exception as e:
            raise RuntimeError(
                f"BLOCKER: TreeSBM project imports failed: {e}"
            ) from e

        try:
            self._torch = torch
            self._mutate = mutate_sequence_independent
            self._get_lm_logits = get_lm_logits
            self._tilt_log_R0 = tilt_log_R0_by_fitness
            self._apply_shm = apply_shm_site_prior
            self._build_seq_indices = _build_seq_indices
            self._TreeState = TreeState
            self._build_edges = build_edges
            self._compute_laplacian_pe = compute_laplacian_pe
            self._compute_structural_features = compute_structural_features

            self.node_enc, self.tree_enc, self.rate_heads, self.col_entropy = load_models(
                self.checkpoint, self.device, self.max_seq_len
            )
            self.embedder = ESM2Embedder(device=self.device)
            # Q0: ESM-2-8M for pathogen / OAS v1/v2. Thrifty AA Q0 only when the
            # ckpt was trained with --r0-backend thrifty_aa (Ab Recipe B).
            cfg = getattr(self.rate_heads, "_ckpt_config", None) or {}
            r0_name = normalize_backend_name(cfg.get("r0_backend") or "esm2")
            self._r0_name = r0_name
            self._r0 = None
            self.tokenizer = None
            self.esm_model = None
            self.aa_token_ids = None
            if r0_name == BACKEND_THRIFTY_AA:
                self._r0 = build_r0_backend(r0_name, device=self.device)
            else:
                esm_id = "facebook/esm2_t6_8M_UR50D"
                self.tokenizer = AutoTokenizer.from_pretrained(esm_id)
                self.esm_model = EsmForMaskedLM.from_pretrained(esm_id).to(self.device).eval()
                for p in self.esm_model.parameters():
                    p.requires_grad = False
                self.aa_token_ids = torch.tensor(
                    [self.tokenizer.convert_tokens_to_ids(aa) for aa in AA_VOCAB],
                    dtype=torch.long,
                )
            self._ready = True
        except Exception as e:
            raise RuntimeError(
                f"BLOCKER: failed to load TreeSBM checkpoint {ckpt}: {e}"
            ) from e

    def sample_child(self, parent_seq: str, branch_length: float, rng_seed: int) -> str:
        if not self._ready:
            raise RuntimeError("TreeSBM sampler not loaded")
        torch = self._torch
        torch.manual_seed(int(rng_seed))

        # Minimal single-node tree: evolve the parent sequence for observed BL.
        tree = self._TreeState.root_only(parent_seq)
        node_ids_t = tree.node_ids
        node_to_idx = {nid: i for i, nid in enumerate(node_ids_t)}
        active_leaves = list(tree.active_leaves)
        active_idx = [node_to_idx[v] for v in active_leaves]
        node_times_dict = {nid: 0.0 for nid in node_ids_t}

        struct_t = self._compute_structural_features(tree, node_to_idx).to(self.device)
        lap_t = self._compute_laplacian_pe(tree, node_to_idx, 8, device=self.device)
        edge_index_t, _, edge_attr_t = self._build_edges(tree, node_to_idx)
        edge_index_t = edge_index_t.to(self.device)
        branch_lens_t = edge_attr_t.squeeze(-1).to(self.device)
        node_seqs = [tree.node_seqs[nid] for nid in node_ids_t]
        plm_t = self.embedder.embed_sequences(node_seqs).to(self.device)
        active_seqs = [tree.node_seqs[v] for v in active_leaves]
        if self._r0 is not None:
            log_R0_mut = self._r0.log_mutation_rates(
                active_seqs, self.max_seq_len, device=self.device
            )
            if not torch.is_tensor(log_R0_mut):
                raise TypeError("r0_backend.log_mutation_rates must return a tensor")
            log_R0_mut = log_R0_mut.to(self.device)
        else:
            log_R0_mut = self._get_lm_logits(
                self.tokenizer,
                self.esm_model,
                self.aa_token_ids,
                active_seqs,
                self.max_seq_len,
                self.device,
            )
        log_R0_mut = self._tilt_log_R0(
            log_R0_mut,
            beta=float(getattr(self.rate_heads, "_fitness_beta", 0.0) or 0.0),
            score=getattr(self.rate_heads, "_fitness_score", "log_R0"),
            mode=getattr(self.rate_heads, "_fitness_tilt_mode", "site_local"),
        )
        shm_boost = float(getattr(self.rate_heads, "_shm_site_boost", 0.0) or 0.0)
        shm_fwr = float(getattr(self.rate_heads, "_shm_fwr_stay", 0.0) or 0.0)
        shm_aid = bool(getattr(self.rate_heads, "_shm_use_aid", False))
        if shm_boost != 0.0 or shm_fwr != 0.0 or shm_aid:
            log_R0_mut = self._apply_shm(
                log_R0_mut,
                active_seqs,
                cdr_mask=getattr(self.rate_heads, "_mut_hotspot_mask", None),
                boost=shm_boost,
                fwr_stay=shm_fwr,
                use_aid=shm_aid,
            )
        log_pssm = getattr(self.rate_heads, "_train_log_pssm", None)
        aa_indices = None
        if getattr(self.rate_heads, "needs_aa_indices", False):
            # OAS / mut-aa-emb / cosine-head ckpts need current-AA indices.
            aa_indices = self._build_seq_indices(
                active_seqs, self.max_seq_len, self.device
            )
        with torch.no_grad():
            h_t = self.node_enc(plm_t, struct_t, lap_t)
            H_t, _ = self.tree_enc(
                h_t,
                node_ids_t,
                node_times_dict,
                edge_index_t,
                branch_lens_t,
                t_scalar=0.0,
            )
            out = self.rate_heads(
                H_t,
                active_idx,
                log_R0_mut,
                site_entropy=self.col_entropy,
                aa_indices=aa_indices,
                log_pssm=log_pssm,
            )
        log_R_i = out["log_R_theta_mut"][0]
        seq_len = min(len(parent_seq), self.max_seq_len)
        return self._mutate(
            log_R_i,
            parent_seq,
            seq_len,
            float(branch_length),
            self.mutation_rate_scale,
        )
