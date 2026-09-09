"""Appendix E.2 mut-head / stop-head train flags."""

from pathlib import Path

import torch

from src.networks import RateHeads

ROOT = Path(__file__).resolve().parents[1]


def test_train_py_has_e2_flags():
    src = (ROOT / "scripts/train.py").read_text()
    assert "--ablate-mut-head" in src
    assert "--ablate-stop-head" in src
    assert "ablate_mut_head" in src
    assert "ablate_stop_head" in src


def test_load_models_reads_e2_ckpt_config():
    src = (ROOT / "scripts/eval_single_tree.py").read_text()
    assert 'cfg.get("ablate_mut_head"' in src
    assert 'cfg.get("ablate_stop_head"' in src


def test_ablate_mut_head_returns_r0():
    heads = RateHeads(d_model=128, max_seq_len=32)
    heads.ablate_mut_head = True
    H = torch.randn(3, 128)
    log_R0 = torch.randn(2, 32, 20)
    out = heads(H, [0, 1], log_R0)
    assert torch.allclose(out["log_R_theta_mut"], log_R0)
    assert out["branching_rate"].shape == (2,)
    assert out["branch_length"].shape == (2,)


def test_ablate_stop_head_constant():
    heads = RateHeads(d_model=128, max_seq_len=32)
    heads.ablate_stop_head = True
    H = torch.randn(3, 128)
    log_R0 = torch.randn(2, 32, 20)
    out = heads(H, [0, 1], log_R0)
    assert torch.allclose(out["stop_prob"], torch.full((2,), 0.5))
    assert out["log_R_theta_mut"].shape == (2, 32, 20)
