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


def select_mut_hotspots(
    site_entropy: torch.Tensor,
    topk: int | None = None,
    frac: float | None = None,
) -> torch.Tensor:
    """
    Binary hotspot mask over alignment columns from TRAIN MSA column entropy.

    Selects the N highest-entropy sites (``topk``) or the top ``frac`` fraction
    of columns (``ceil(frac * L)``, at least 1 if frac>0). Exactly one of
    ``topk`` / ``frac`` must be set. Ties are broken by column index (stable).

    Returns a bool tensor of shape [L]. Off-by-default callers pass this into
    ``bridge_losses`` to hard-boost L_mut at mutating regions (vs soft
    floor+alpha*H, which reweights all columns continuously).
    """
    entropy = torch.as_tensor(site_entropy)
    if entropy.ndim != 1:
        raise ValueError(f"site_entropy for hotspot selection must be [L], got {tuple(entropy.shape)}")
    L = int(entropy.numel())
    if L == 0:
        return torch.zeros(0, dtype=torch.bool, device=entropy.device)

    if (topk is None) == (frac is None):
        raise ValueError("Provide exactly one of topk or frac for hotspot selection")
    if topk is not None:
        if topk < 0:
            raise ValueError("mut_hotspot_topk must be non-negative")
        k = min(int(topk), L)
    else:
        if not (0.0 <= frac <= 1.0):
            raise ValueError("mut_hotspot_frac must be in [0, 1]")
        k = int(math.ceil(frac * L)) if frac > 0 else 0
        k = min(k, L)

    mask = torch.zeros(L, dtype=torch.bool, device=entropy.device)
    if k == 0:
        return mask
    # Stable top-k: highest entropy first; equal entropy → lower index first.
    order = torch.argsort(entropy, descending=True, stable=True)
    mask[order[:k]] = True
    return mask


def _expand_hotspot_mask(
    mut_hotspot_mask: torch.Tensor,
    n: int,
    max_seq_len: int,
    device,
) -> torch.Tensor:
    """Broadcast hotspot mask to [n, L] bool."""
    hot = torch.as_tensor(mut_hotspot_mask, device=device)
    if hot.ndim == 1:
        hot = hot.unsqueeze(0)
    if hot.ndim != 2:
        raise ValueError(
            "mut_hotspot_mask must have shape [L], [1, L], or [n, L], "
            f"got {tuple(hot.shape)}"
        )
    if hot.shape[1] != max_seq_len:
        raise ValueError(
            f"mut_hotspot_mask length {hot.shape[1]} does not match "
            f"max_seq_len={max_seq_len}"
        )
    if hot.shape[0] == 1:
        hot = hot.expand(n, -1)
    elif hot.shape[0] != n:
        raise ValueError(
            f"mut_hotspot_mask batch dimension must be 1 or {n}, got {hot.shape[0]}"
        )
    return hot.bool()


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
    lambda_cons: float = 1.0,
    bridge_c: float = 1.0,
    device: str = "cpu",
    site_entropy: torch.Tensor | None = None,
    use_entropy_loss_weighting: bool = False,
    use_entropy_cons_weighting: bool = False,
    entropy_weight_alpha: float = 1.0,
    entropy_weight_alpha_cons: float | None = None,
    entropy_weight_floor: float = 1.0,
    entropy_is_normalized: bool = False,
    mut_normalize: str = "mean",
    mut_hotspot_mask: torch.Tensor | None = None,
    mut_hotspot_weight: float = 1.0,
    mut_hotspot_force: bool = False,
) -> dict:
    """
    Bridge matching losses.

    Hard MSA-entropy hotspots (``mut_hotspot_mask`` from ``select_mut_hotspots``):
      Default: multiply L_mut site weights by ``mut_hotspot_weight`` on
      mut_mask ∩ hotspot (soft floor+αH still applies if enabled). Does **not**
      pull conserved (aa_t==x1) positions into L_mut.
      Optional ``mut_hotspot_force``: L_mut mask becomes
      (mut_mask | hotspot) & valid, and those hotspot∩cons sites are removed from
      cons_mask so high-entropy columns train as mutation sites even when this
      sample is already at the T1 AA.
    """
    n = len(active_leaves)
    eps_rate = 1e-6

    if entropy_weight_alpha < 0:
        raise ValueError("entropy_weight_alpha must be non-negative")
    if entropy_weight_alpha_cons is not None and entropy_weight_alpha_cons < 0:
        raise ValueError("entropy_weight_alpha_cons must be non-negative")
    if entropy_weight_floor <= 0:
        raise ValueError("entropy_weight_floor must be greater than zero")
    if lambda_cons < 0:
        raise ValueError("lambda_cons must be non-negative")
    if mut_normalize not in ("mean", "count"):
        raise ValueError("mut_normalize must be 'mean' or 'count'")
    if mut_hotspot_weight <= 0:
        raise ValueError("mut_hotspot_weight must be greater than zero")
    if mut_hotspot_force and mut_hotspot_mask is None:
        raise ValueError("mut_hotspot_force requires mut_hotspot_mask")

    alpha_cons = (
        entropy_weight_alpha
        if entropy_weight_alpha_cons is None
        else entropy_weight_alpha_cons
    )

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
            "L_semi": z,
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
    # Soft alignment-entropy weighting (optional):
    #   L_mut  weight = floor + alpha * entropy         -> mutate freely at hotspots
    #   L_cons weight = floor + alpha * (1 - entropy)    -> stay put at cold sites
    # Hard top-k/frac hotspot boost (optional, stacks with soft):
    #   L_mut weight *= mut_hotspot_weight on hotspot columns (see docstring).
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

    hotspot_2d = None
    if mut_hotspot_mask is not None:
        hotspot_2d = _expand_hotspot_mask(
            mut_hotspot_mask, n, max_seq_len, kl_per_pos.device
        )
        if mut_hotspot_force:
            # High-entropy columns contribute to L_mut even when aa_t == x1
            # this step; drop them from cons_mask to avoid opposing gradients.
            forced = hotspot_2d & valid_mask
            mut_mask = mut_mask | forced
            cons_mask = cons_mask & ~forced

    n_mut = mut_mask.sum().clamp_min(1).to(dtype=kl_per_pos.dtype)

    # Force only expands mut_mask; weight≠1 enters the weighted path.
    use_weighted_mut = use_entropy_loss_weighting or (
        hotspot_2d is not None and mut_hotspot_weight != 1.0
    )
    if use_weighted_mut:
        if use_entropy_loss_weighting:
            site_weights = entropy_weight_floor + entropy_weight_alpha * normalized_entropy
        else:
            site_weights = torch.ones_like(kl_per_pos)
        if hotspot_2d is not None and mut_hotspot_weight != 1.0:
            site_weights = torch.where(
                hotspot_2d,
                site_weights * mut_hotspot_weight,
                site_weights,
            )
        if mut_mask.any():
            mut_kl = kl_per_pos[mut_mask]
            mut_weights = site_weights[mut_mask]
            weighted_sum = (mut_kl * mut_weights).sum()
            # mean: antiGen-style /Z with Z = sum(weights); count: /n_mut so hotspot
            # weights increase total mut mass rather than only rebalancing within muts.
            if mut_normalize == "count":
                L_mut = weighted_sum / n_mut
            else:
                L_mut = weighted_sum / mut_weights.sum().clamp_min(1e-8)
            if normalized_entropy is not None:
                mean_mut_entropy = normalized_entropy[mut_mask].detach().mean()
            else:
                mean_mut_entropy = torch.zeros((), device=kl_per_pos.device)
            mean_mut_weight = mut_weights.detach().mean()
            max_mut_weight = mut_weights.detach().max()
        else:
            L_mut = kl_per_pos.sum() * 0.0
            mean_mut_entropy = torch.zeros((), device=kl_per_pos.device)
            mean_mut_weight = torch.zeros((), device=kl_per_pos.device)
            max_mut_weight = torch.zeros((), device=kl_per_pos.device)
    else:
        if mut_mask.any():
            # mean and count coincide without entropy weights (both /n_mut).
            L_mut = kl_per_pos[mut_mask].sum() / n_mut
        else:
            L_mut = kl_per_pos.sum() * 0.0
        mean_mut_entropy = torch.zeros((), device=kl_per_pos.device)
        mean_mut_weight = torch.ones((), device=kl_per_pos.device)
        max_mut_weight = torch.ones((), device=kl_per_pos.device)

    if use_entropy_cons_weighting and cons_mask.any():
        cons_weights = entropy_weight_floor + alpha_cons * (1.0 - normalized_entropy)
        cons_kl = kl_per_pos[cons_mask]
        cw = cons_weights[cons_mask]
        L_cons = (cons_kl * cw).sum() / cw.sum().clamp_min(1e-8)
    else:
        L_cons = kl_per_pos[cons_mask].mean() if cons_mask.any() else kl_per_pos.sum() * 0.0
    L_rate = lambda_mut * L_mut + lambda_cons * L_cons

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

    # L_semi is optional; callers add it via attach_semigroup_loss when λ_semi > 0.
    L_semi = torch.zeros((), device=device)
    total = (
        L_rate
        + lambda_top * L_top
        + lambda_br * L_br
        + lambda_stop * L_stop
        + lambda_pll * L_pll
    )
    return {
        "L_rate": L_rate, "L_mut": L_mut, "L_cons": L_cons,
        "L_top": L_top, "L_br": L_br, "L_stop": L_stop, "L_pll": L_pll,
        "L_semi": L_semi,
        "L_br_pred_std": br_pred_std, "L_br_target_std": br_target_std,
        "mean_mut_entropy": mean_mut_entropy,
        "mean_mut_weight": mean_mut_weight,
        "max_mut_weight": max_mut_weight,
        "total": total,
    }


def attach_semigroup_loss(
    losses: dict,
    L_semi: torch.Tensor,
    lambda_semi: float,
) -> dict:
    """Fold λ_semi * L_semi into an existing bridge_losses() result dict."""
    losses = dict(losses)
    losses["L_semi"] = L_semi
    if lambda_semi != 0.0:
        losses["total"] = losses["total"] + lambda_semi * L_semi
    return losses
