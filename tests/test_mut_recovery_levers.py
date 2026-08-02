"""Unit tests for mut_recovery train/inference levers (SH26 iterations)."""

import torch

from src.bridge.losses import bridge_losses
from src.bridge.mutation_sample import (
    AA_VOCAB,
    mutate_sequence_independent,
    mutate_sequence_site_softmax,
)
from src.networks import RateHeads


def _toy_bridge_inputs(n=2, L=4, device="cpu"):
    log_R = torch.randn(n, L, 20, device=device)
    log_R0 = torch.randn(n, L, 20, device=device)
    seqs_t = ["ACDE", "ACDE"]
    targets = {"n0": "AADE", "n1": "ACDF"}  # 1 mut each
    active = ["n0", "n1"]
    child_counts = {nid: 0 for nid in active}
    child_bls = {nid: [] for nid in active}
    return dict(
        log_R_theta_mut=log_R,
        log_R_theta_branch=torch.ones(n, device=device),
        branch_length_pred=torch.ones(n, device=device) * 0.1,
        stop_prob=torch.ones(n, device=device) * 0.5,
        log_R0_mut=log_R0,
        seqs_t=seqs_t,
        active_leaves=active,
        T1_mut_targets=targets,
        T1_child_counts=child_counts,
        T1_child_bls=child_bls,
        t=0.5,
        max_seq_len=L,
        device=device,
    )


def test_lambda_cons_scales_l_rate():
    kwargs = _toy_bridge_inputs()
    base = bridge_losses(**kwargs, lambda_mut=1.0, lambda_cons=1.0)
    down = bridge_losses(**kwargs, lambda_mut=1.0, lambda_cons=0.25)
    expected = base["L_mut"] + 0.25 * base["L_cons"]
    assert torch.allclose(down["L_rate"], expected, atol=1e-5)


def test_mut_normalize_count_vs_mean_with_entropy():
    kwargs = _toy_bridge_inputs(L=4)
    site_entropy = torch.tensor([0.1, 0.9, 0.1, 0.1])
    mean_out = bridge_losses(
        **kwargs,
        site_entropy=site_entropy,
        use_entropy_loss_weighting=True,
        entropy_weight_alpha=3.0,
        entropy_is_normalized=True,
        mut_normalize="mean",
    )
    count_out = bridge_losses(
        **kwargs,
        site_entropy=site_entropy,
        use_entropy_loss_weighting=True,
        entropy_weight_alpha=3.0,
        entropy_is_normalized=True,
        mut_normalize="count",
    )
    assert mean_out["L_mut"].ndim == 0
    assert count_out["L_mut"].ndim == 0
    assert not torch.allclose(mean_out["L_mut"], count_out["L_mut"], atol=1e-6)


def test_entropy_alpha_cons_independent():
    kwargs = _toy_bridge_inputs()
    site_entropy = torch.tensor([0.0, 1.0, 0.0, 0.0])
    hi = bridge_losses(
        **kwargs,
        site_entropy=site_entropy,
        use_entropy_cons_weighting=True,
        entropy_weight_alpha=1.0,
        entropy_weight_alpha_cons=5.0,
        entropy_is_normalized=True,
    )
    lo = bridge_losses(
        **kwargs,
        site_entropy=site_entropy,
        use_entropy_cons_weighting=True,
        entropy_weight_alpha=1.0,
        entropy_weight_alpha_cons=0.0,
        entropy_is_normalized=True,
    )
    assert not torch.allclose(hi["L_cons"], lo["L_cons"], atol=1e-6)


def test_rateheads_default_matches_legacy_mut_in():
    """Flag-off RateHeads keeps mut_in = d_model+20 (checkpoint-compatible)."""
    heads = RateHeads(d_model=128, max_seq_len=32)
    assert not heads.use_mut_aa_emb
    assert not heads.use_pssm_gate
    H = torch.randn(3, 128)
    log_R0 = torch.randn(2, 32, 20)
    out = heads(H, [0, 2], log_R0)
    assert out["log_R_theta_mut"].shape == (2, 32, 20)


def test_rateheads_mut_aa_emb_and_pssm_gate():
    L = 16
    heads = RateHeads(
        d_model=64, max_seq_len=L,
        use_mut_aa_emb=True, d_aa=8,
        use_pssm_gate=True, pssm_gate_fixed_w=0.75,
    )
    H = torch.randn(2, 64)
    log_R0 = torch.randn(1, L, 20)
    aa = torch.randint(0, 20, (1, L))
    pssm = torch.randn(L, 20)
    out = heads(H, [0], log_R0, aa_indices=aa, log_pssm=pssm)
    assert out["log_R_theta_mut"].shape == (1, L, 20)
    assert torch.isfinite(out["log_R_theta_mut"]).all()


def test_legacy_state_dict_loads_without_new_flags():
    legacy = RateHeads(d_model=64, max_seq_len=16)
    sd = legacy.state_dict()
    fresh = RateHeads(d_model=64, max_seq_len=16)
    fresh.load_state_dict(sd)


def test_site_softmax_and_independent_samplers():
    torch.manual_seed(0)
    L = 8
    log_R = torch.zeros(L, 20)
    w_idx = AA_VOCAB.index("W")
    log_R[3, :] = -10.0
    log_R[3, w_idx] = 5.0
    seq = "A" * L
    out = mutate_sequence_site_softmax(log_R, seq, L, dt=1.0, mutation_rate_scale=5.0)
    out2 = mutate_sequence_independent(log_R, seq, L, dt=1.0, mutation_rate_scale=5.0)
    assert len(out) == L
    assert len(out2) == L
