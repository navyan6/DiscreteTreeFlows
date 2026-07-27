"""
Bridge matching training losses

L_rate: Algorithm 1 bridge matching — KL( R^{0|T1}_t || R_theta ) on the mutation
        head. The target is the reference process P^0 Doob h-transformed to hit the
        observed terminal AA x1 (see src/bridge/conditional_rates.py). Recovers the
        old terminal cross-entropy as the t->1 limit; for t<1 it anchors off-target
        mass to the ESM reference q = softmax(log_R0).
L_top:  Poisson NLL on branching rate vs. T1 child count
importnant note: fasttree assumes a bifurcating tree, so each parent can either have 0, 1, or max two chldren
L_br:   MSE on predicted branch length vs. mean T1 child branch length
L_stop: BCE on stop_prob vs. whether leaf has no children in T1
L_pll:  ESM PLL regularizer - penalizes sequences drifting from ESM fitness landscape (we dont want nonsensical sequences!)
total:  weighted sum of all five terms
"""

import math

import torch
import torch.nn.functional as F

from src.bridge.conditional_rates import conditional_bridge_kl

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_VOCAB)}
PAD_IDX = 20


def _build_aa_targets(active_leaves, T1_mut_targets, max_seq_len, device):
    # Leaves absent from T1_mut_targets (terminal T1 leaves) stay all-PAD → excluded by valid_mask
    n = len(active_leaves)
    targets = torch.full((n, max_seq_len), PAD_IDX, dtype=torch.long, device=device)
    for i, nid in enumerate(active_leaves):
        seq = T1_mut_targets.get(nid, "")
        for j, aa in enumerate(seq[:max_seq_len]):
            targets[i, j] = AA_TO_IDX.get(aa, PAD_IDX)
    return targets


def _build_seq_indices(seqs_t, max_seq_len, device):
    #Convert current T_t sequences to integer index tensor [n, max_seq_len].
    n = len(seqs_t)
    indices = torch.full((n, max_seq_len), PAD_IDX, dtype=torch.long, device=device)
    for i, seq in enumerate(seqs_t):
        for j, aa in enumerate(seq[:max_seq_len]):
            indices[i, j] = AA_TO_IDX.get(aa, PAD_IDX)
    return indices


def _prepare_alignment_entropy(
    site_entropy: torch.Tensor,
    n: int,
    max_seq_len: int,
    device,
    dtype: torch.dtype,
    entropy_is_normalized: bool,
) -> torch.Tensor:
    """
    Prepare per-site Shannon entropy computed from the TRAINING alignment only.

    Accepted shapes:
      - [L]
      - [1, L]
      - [n, L]
    """
    entropy = torch.as_tensor(site_entropy, dtype=dtype, device=device)

    if entropy.ndim == 1:
        entropy = entropy.unsqueeze(0)

    if entropy.ndim != 2:
        raise ValueError(
            "site_entropy must have shape [L], [1, L], or [n, L], "
            f"got {tuple(entropy.shape)}"
        )

    if entropy.shape[1] != max_seq_len:
        raise ValueError(
            f"site_entropy length {entropy.shape[1]} does not match "
            f"max_seq_len={max_seq_len}"
        )

    if entropy.shape[0] == 1:
        entropy = entropy.expand(n, -1)
    elif entropy.shape[0] != n:
        raise ValueError(
            f"site_entropy batch dimension must be 1 or {n}, got {entropy.shape[0]}"
        )

    if not torch.isfinite(entropy).all():
        raise ValueError("site_entropy contains NaN or infinite values")
    if (entropy < 0).any():
        raise ValueError("site_entropy must be non-negative")

    if not entropy_is_normalized:
        entropy = entropy / math.log(len(AA_VOCAB))

    return entropy.clamp(min=0.0, max=1.0)


def bridge_losses(
    log_R_theta_mut: torch.Tensor,
    log_R_theta_branch: torch.Tensor,
    branch_length_pred: torch.Tensor,
    stop_prob: torch.Tensor,
    log_R0_mut: torch.Tensor | None,
    seqs_t: list[str],
    active_leaves: list[str],
    T1_mut_targets: dict[str, str],
    T1_child_counts: dict[str, int],
    T1_child_bls: dict[str, list[float]],
    t: float,
    max_seq_len: int,
    lambda_top: float = 0.1,
    lambda_br: float = 0.1,
    lambda_stop: float = 0.1,
    lambda_pll: float = 0.01,
    lambda_mut: float = 5.0,
    bridge_c: float = 1.0,
    device: str = "cpu",
    site_entropy: torch.Tensor | None = None,
    use_entropy_loss_weighting: bool = False,
    use_entropy_cons_weighting: bool = False,
    entropy_weight_alpha: float = 1.0,
    entropy_weight_floor: float = 1.0,
    entropy_is_normalized: bool = False,
) -> dict:
    n = len(active_leaves)
    eps_rate = 1e-6

    if entropy_weight_alpha < 0:
        raise ValueError("entropy_weight_alpha must be non-negative")
    if entropy_weight_floor <= 0:
        raise ValueError("entropy_weight_floor must be greater than zero")

    if n == 0:
        z = torch.zeros((), device=device, requires_grad=True)
        return {
            "L_rate": z,
            "L_mut": z,
            "L_cons": z,
            "L_top": z,
            "L_br": z,
            "L_stop": z,
            "L_pll": z,
            "L_br_pred_std": z.detach(),
            "L_br_target_std": z.detach(),
            "mean_mut_entropy": z.detach(),
            "mean_mut_weight": z.detach(),
            "max_mut_weight": z.detach(),
            "total": z,
        }

    # ── L_rate: Algorithm-1 bridge matching, KL( R^{0|T1}_t || R_theta )
    # Target = reference P^0 Doob h-transformed to the terminal AA x1 (conditional_rates).
    targets  = _build_aa_targets(active_leaves, T1_mut_targets, max_seq_len, device)  # [n, L] x1 (sampled T1 leaf AAs)
    aa_t     = _build_seq_indices(seqs_t, max_seq_len, device)                  # [n, L] a  (T_t AAs)

    ref_logits = log_R0_mut if log_R0_mut is not None else torch.zeros_like(log_R_theta_mut)
    kl_per_pos = conditional_bridge_kl(
        log_R_theta_mut, ref_logits, targets, t=t, c=bridge_c
    )                                                                           # [n, L]

    valid_mask = (targets != PAD_IDX) & (aa_t != PAD_IDX)
    mut_mask   = (aa_t != targets) & valid_mask   # positions that mutate T_t→T1
    cons_mask  = (aa_t == targets) & valid_mask   # positions already at T1 AA

    # Upweight rare mutating positions (sparse signal); time-weighting is already
    # handled inside the h-transform, so no extra 1/(1-t) factor.
    #
    # Optional alignment-entropy weighting (computed once, used by both terms):
    #   L_mut  weight = floor + alpha * entropy         -> mutate freely at hotspots
    #   L_cons weight = floor + alpha * (1 - entropy)    -> stay put at cold sites
    # The L_cons variant directly fights over-mutation of conserved regions.
    normalized_entropy = None
    if use_entropy_loss_weighting or use_entropy_cons_weighting:
        if site_entropy is None:
            raise ValueError(
                "site_entropy is required when use_entropy_loss_weighting or "
                "use_entropy_cons_weighting is True."
            )
        normalized_entropy = _prepare_alignment_entropy(
            site_entropy=site_entropy,
            n=n,
            max_seq_len=max_seq_len,
            device=kl_per_pos.device,
            dtype=kl_per_pos.dtype,
            entropy_is_normalized=entropy_is_normalized,
        )

    if use_entropy_loss_weighting:
        site_weights = entropy_weight_floor + entropy_weight_alpha * normalized_entropy
        if mut_mask.any():
            mut_kl = kl_per_pos[mut_mask]
            mut_weights = site_weights[mut_mask]
            L_mut = (mut_kl * mut_weights).sum() / mut_weights.sum().clamp_min(1e-8)
            mean_mut_entropy = normalized_entropy[mut_mask].detach().mean()
            mean_mut_weight = mut_weights.detach().mean()
            max_mut_weight = mut_weights.detach().max()
        else:
            L_mut = kl_per_pos.sum() * 0.0
            mean_mut_entropy = torch.zeros((), device=kl_per_pos.device)
            mean_mut_weight = torch.zeros((), device=kl_per_pos.device)
            max_mut_weight = torch.zeros((), device=kl_per_pos.device)
    else:
        L_mut = kl_per_pos[mut_mask].mean() if mut_mask.any() else kl_per_pos.sum() * 0.0
        mean_mut_entropy = torch.zeros((), device=kl_per_pos.device)
        mean_mut_weight = torch.ones((), device=kl_per_pos.device)
        max_mut_weight = torch.ones((), device=kl_per_pos.device)

    if use_entropy_cons_weighting and cons_mask.any():
        cons_weights = entropy_weight_floor + entropy_weight_alpha * (1.0 - normalized_entropy)
        cons_kl = kl_per_pos[cons_mask]
        cw = cons_weights[cons_mask]
        L_cons = (cons_kl * cw).sum() / cw.sum().clamp_min(1e-8)
    else:
        L_cons = kl_per_pos[cons_mask].mean() if cons_mask.any() else kl_per_pos.sum() * 0.0
    L_rate = lambda_mut * L_mut + L_cons

    # ── L_top 
    child_counts = torch.tensor(
        [T1_child_counts[nid] for nid in active_leaves],
        dtype=torch.float32, device=device,
    )
    L_top = F.poisson_nll_loss(
        torch.log(log_R_theta_branch + eps_rate),
        child_counts,
        log_input=True,
        full=False,
    )

    # ── L_br  (plain MSE; log-space + terminal-masking variant was tried and made
    # generated branches worse — 9.5x vs 4.8x — so reverted. Branch-length skew still
    # unsolved; needs per-child log-branch modeling, a separate redesign.)
    target_bls = torch.tensor(
        [
            (sum(T1_child_bls[nid]) / len(T1_child_bls[nid]))
            if T1_child_bls[nid] else 0.0
            for nid in active_leaves
        ],
        dtype=torch.float32, device=device,
    )
    L_br = F.mse_loss(branch_length_pred, target_bls)
    # diagnostics only (not in `total`): tell apart "targets are just tiny" from
    # "prediction collapsed to a near-constant regardless of input"
    br_pred_std = branch_length_pred.detach().std() if branch_length_pred.numel() > 1 \
        else torch.zeros((), device=device)
    br_target_std = target_bls.std() if target_bls.numel() > 1 else torch.zeros((), device=device)

    # ── L_stop 
    has_no_children = torch.tensor(
        [T1_child_counts[nid] == 0 for nid in active_leaves],
        dtype=torch.float32, device=device,
    )
    L_stop = F.binary_cross_entropy(stop_prob, has_no_children)

    # ── L_pll 
    if log_R0_mut is not None:
        aa_indices = _build_seq_indices(seqs_t, max_seq_len, device)
        pll_mask   = aa_indices != PAD_IDX
        # clamp before gather so PAD_IDX=20 doesn't go out-of-bounds on dim size 20
        aa_safe    = aa_indices.clamp(0, 19)
        pll_scores = log_R0_mut.gather(-1, aa_safe.unsqueeze(-1)).squeeze(-1)
        L_pll = -pll_scores[pll_mask].mean()
    else:
        L_pll = torch.zeros((), device=device)

    total = L_rate + lambda_top * L_top + lambda_br * L_br + lambda_stop * L_stop + lambda_pll * L_pll
    return {
        "L_rate": L_rate, "L_mut": L_mut, "L_cons": L_cons,
        "L_top": L_top, "L_br": L_br, "L_stop": L_stop, "L_pll": L_pll,
        "L_br_pred_std": br_pred_std, "L_br_target_std": br_target_std,
        "mean_mut_entropy": mean_mut_entropy,
        "mean_mut_weight": mean_mut_weight,
        "max_mut_weight": max_mut_weight,
        "total": total,
    }
