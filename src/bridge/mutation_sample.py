"""Inference-time AA mutation samplers for controlled tree generation."""

from __future__ import annotations

import torch

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_VOCAB)}


def mutate_sequence_independent(
    log_R_theta: torch.Tensor,
    seq: str,
    seq_len: int,
    dt: float,
    mutation_rate_scale: float = 1.0,
) -> str:
    """
    Per-site independent Bernoulli mutations (default TreeSBM generate path).

    log_R_theta: [L, 20]
    """
    new_seq = list(seq)
    for pos in range(seq_len):
        curr_idx = AA_TO_IDX.get(seq[pos], -1)
        if curr_idx < 0:
            continue
        probs = log_R_theta[pos].softmax(-1)
        probs_mut = probs.clone()
        probs_mut[curr_idx] = 0.0
        total = probs_mut.sum().item()
        if total > 0 and torch.rand(1).item() < total * dt * mutation_rate_scale:
            new_seq[pos] = AA_VOCAB[torch.multinomial(probs_mut / total, 1).item()]
    return "".join(new_seq)


def mutate_sequence_site_softmax(
    log_R_theta: torch.Tensor,
    seq: str,
    seq_len: int,
    dt: float,
    mutation_rate_scale: float = 1.0,
    site_temperature: float = 1.0,
) -> str:
    """
    Site-propensity sampling (antiGen-style joint L×A mass analogue).

    1. Site fire score ∝ Σ_{a≠curr} softmax(log R)_a  (= 1 - p_stay)
    2. n_events ~ Poisson(Σ scores · dt · mrs), clamped to seq_len
    3. For each event: sample site ~ Categorical(scores), then AA | site
    """
    if site_temperature <= 0:
        raise ValueError("site_temperature must be positive")

    new_seq = list(seq)
    site_scores = []
    mut_dists = []
    valid_pos = []

    for pos in range(seq_len):
        curr_idx = AA_TO_IDX.get(seq[pos], -1)
        if curr_idx < 0:
            continue
        probs = log_R_theta[pos].softmax(-1)
        probs_mut = probs.clone()
        probs_mut[curr_idx] = 0.0
        total = probs_mut.sum()
        if total.item() <= 0:
            continue
        valid_pos.append(pos)
        site_scores.append(total)
        mut_dists.append(probs_mut / total)

    if not valid_pos:
        return "".join(new_seq)

    scores = torch.stack(site_scores)  # [S]
    mean_events = scores.sum().item() * dt * mutation_rate_scale
    n_events = int(torch.poisson(torch.tensor(mean_events)).item())
    n_events = max(0, min(n_events, seq_len))
    if n_events == 0:
        return "".join(new_seq)

    logits = (scores.clamp_min(1e-12).log() / site_temperature)
    site_probs = torch.softmax(logits, dim=0)
    for _ in range(n_events):
        j = int(torch.multinomial(site_probs, 1).item())
        pos = valid_pos[j]
        # Recompute vs current AA in case the same site is hit twice in one step.
        curr_idx = AA_TO_IDX.get(new_seq[pos], -1)
        if curr_idx < 0:
            continue
        probs = log_R_theta[pos].softmax(-1)
        probs_mut = probs.clone()
        probs_mut[curr_idx] = 0.0
        total = probs_mut.sum()
        if total.item() <= 0:
            continue
        new_seq[pos] = AA_VOCAB[torch.multinomial(probs_mut / total, 1).item()]

    return "".join(new_seq)
