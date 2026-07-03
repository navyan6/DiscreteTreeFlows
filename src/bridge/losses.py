"""
Bridge matching training losses (Algorithm 1 of TreeSBM).

L_seq: cross-entropy on predicting T1 AA sequence at active leaves,
       time-weighted by 1/(1-t) (Doob h-transform)
L_top: Poisson NLL on branching rate vs. actual T1 child count
L_br:  MSE on predicted branch length vs. mean T1 child branch length
total: L_seq + lambda_top * L_top + lambda_br * L_br
"""

import torch
import torch.nn.functional as F

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_VOCAB)}
PAD_IDX = 20


def bridge_losses(
    out: dict,
    active_leaves: list[str],
    T1_seqs: dict[str, str],
    T1_child_counts: dict[str, int],
    T1_child_bls: dict[str, list[float]],
    t: float,
    max_seq_len: int,
    lambda_top: float = 0.1,
    lambda_br: float = 0.1,
    device: str = "cpu",
) -> dict:
    """
    Args:
        out:               RateHeads output dict (mutation_logits, branching_rate, branch_length)
        active_leaves:     list of active leaf node IDs in T_t (ordered to match out tensors)
        T1_seqs:           T1 AA sequences (prediction targets)
        T1_child_counts:   n_children in T1 per active leaf
        T1_child_bls:      branch lengths to T1 children per active leaf
        t:                 current interpolation time
        max_seq_len:       sequence length (same as RateHeads.max_seq_len)
        lambda_top/br:     auxiliary loss weights
        device:            torch device string
    """
    n = len(active_leaves)
    eps_t = 1e-2
    eps_rate = 1e-6

    if n == 0:
        z = torch.zeros((), device=device, requires_grad=True)
        return {"L_seq": z, "L_top": z, "L_br": z, "total": z}

    # ── L_seq ─────────────────────────────────────────────────────────────────
    logits = out["mutation_logits"]  # [n, max_seq_len, 20]

    targets = torch.full((n, max_seq_len), PAD_IDX, dtype=torch.long, device=device)
    for i, nid in enumerate(active_leaves):
        seq = T1_seqs.get(nid, "")
        for j, aa in enumerate(seq[:max_seq_len]):
            targets[i, j] = AA_TO_IDX.get(aa, PAD_IDX)

    L_seq = F.cross_entropy(
        logits.reshape(n * max_seq_len, 20),
        targets.reshape(n * max_seq_len),
        ignore_index=PAD_IDX,
    )
    # Doob h-transform time-weighting: loss diverges near t=1, clip with eps_t
    L_seq = L_seq / (1.0 - t + eps_t)

    # ── L_top ─────────────────────────────────────────────────────────────────
    br_rate = out["branching_rate"]  # [n], Softplus output > 0
    child_counts = torch.tensor(
        [T1_child_counts[nid] for nid in active_leaves],
        dtype=torch.float32, device=device,
    )
    # Poisson NLL with log-rate input
    L_top = F.poisson_nll_loss(
        torch.log(br_rate + eps_rate),
        child_counts,
        log_input=True,
        full=False,
    )

    # ── L_br ──────────────────────────────────────────────────────────────────
    bl_pred = out["branch_length"]  # [n], Softplus output > 0
    target_bls = torch.tensor(
        [
            (sum(T1_child_bls[nid]) / len(T1_child_bls[nid]))
            if T1_child_bls[nid] else 0.0
            for nid in active_leaves
        ],
        dtype=torch.float32, device=device,
    )
    L_br = F.mse_loss(bl_pred, target_bls)

    total = L_seq + lambda_top * L_top + lambda_br * L_br
    return {"L_seq": L_seq, "L_top": L_top, "L_br": L_br, "total": total}
