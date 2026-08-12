"""Unit tests for multi-pLM / substitution R0 backends + Poisson λ helpers."""

from __future__ import annotations

import numpy as np
import torch

from src.r0_backends import (
    BACKEND_ESM2,
    BACKEND_ESMC,
    BACKEND_PROGEN2,
    STUB_BACKENDS,
    SubstitutionMatrixR0Backend,
    build_r0_backend,
    cache_tag_for_backend,
    list_backends,
    normalize_backend_name,
    ref_rates_filename,
)
from src.reference_process import (
    BranchingIntensityMLP,
    ConstantBranchingIntensity,
    ReferenceProcess,
    sample_poisson_offspring,
)


class FakeR0Backend:
    """Deterministic R0 for ReferenceProcess tests (no GPU / HF)."""

    name = "fake"

    def log_mutation_rates(self, sequences, max_seq_len, device=None):
        N = len(sequences)
        L = max_seq_len
        # Mild preference for AA index 0 (A).
        logits = torch.full((N, L, 20), -2.0)
        logits[:, :, 0] = 1.0
        return torch.log_softmax(logits, dim=-1)

    def close(self):
        return None


def test_normalize_and_cache_tags():
    assert normalize_backend_name("ESM-C") == BACKEND_ESMC
    assert normalize_backend_name("esm2-650m") == "esm2_650m"
    assert cache_tag_for_backend("esm2") == ""
    assert cache_tag_for_backend("esmc") == "_esmc"
    assert cache_tag_for_backend("jtt") == "_jtt"
    assert ref_rates_filename(3, "") == "group_003_ref_rates.pt"
    assert ref_rates_filename(3, "_esmc") == "group_003_ref_rates_esmc.pt"
    assert "esm2" in list_backends()


def test_progen2_stub_raises():
    backend = build_r0_backend(BACKEND_PROGEN2)
    assert BACKEND_PROGEN2 in STUB_BACKENDS
    try:
        backend.log_mutation_rates(["ACDE"], max_seq_len=4)
        assert False, "expected NotImplementedError"
    except NotImplementedError as e:
        assert "stubbed" in str(e).lower() or "ProGen2" in str(e)


def test_neutral_substitution_shape_and_normalized():
    backend = SubstitutionMatrixR0Backend(model="neutral", stay_mass=0.5)
    seqs = ["ACDE", "GGGG"]
    log_R0 = backend.log_mutation_rates(seqs, max_seq_len=6)
    assert log_R0.shape == (2, 6, 20)
    # Valid positions renormalize; gap pad rows stay 0 (untouched).
    for i, seq in enumerate(seqs):
        for pos in range(len(seq)):
            logsum = torch.logsumexp(log_R0[i, pos], dim=-1)
            assert torch.allclose(logsum, torch.tensor(0.0), atol=1e-5)
            # Stay mass on current AA should be highest (0.5).
            aa = seq[pos]
            from src.r0_backends import AA_TO_IDX
            a = AA_TO_IDX[aa]
            probs = log_R0[i, pos].exp()
            assert probs[a].item() == max(probs.tolist())


def test_poisson_offspring_nonneg_and_mean():
    rng = np.random.default_rng(0)
    samples = [sample_poisson_offspring(2.0, 0.5, rng) for _ in range(2000)]
    assert all(s >= 0 for s in samples)
    # E[Poisson(1.0)] ≈ 1
    assert abs(float(np.mean(samples)) - 1.0) < 0.15


def test_reference_rollout_with_fake_r0():
    ref = ReferenceProcess(
        r0_backend=FakeR0Backend(),
        branching_intensity=ConstantBranchingIntensity(2.0),
        beta=0.0,
        rng=np.random.default_rng(1),
    )
    tree = ref.rollout("ACDEFGHIKL", horizon=0.3, dt=0.1, max_nodes=64)
    assert tree.n_nodes() >= 1
    assert tree.root_id in tree.node_seqs
    assert len(tree.node_seqs[tree.root_id]) == 10


def test_fitness_tilt_changes_fake_rollout_rates():
    fake = FakeR0Backend()
    ref0 = ReferenceProcess(fake, ConstantBranchingIntensity(0.01), beta=0.0)
    ref1 = ReferenceProcess(fake, ConstantBranchingIntensity(0.01), beta=2.0)
    p0 = ref0.get_mutation_rates("AAAA")
    p1 = ref1.get_mutation_rates("AAAA")
    # β>0 further concentrates mass on AA 0 (already preferred by FakeR0).
    assert p1[0, 0] > p0[0, 0]


def test_branching_mlp_alias_esm_c_dim():
    mlp = BranchingIntensityMLP(esm_c_dim=16)
    out = mlp(torch.randn(3, 16))
    assert out.shape == (3,)
    assert torch.all(out > 0)


def test_build_unknown_backend_raises():
    try:
        build_r0_backend("not_a_real_backend")
        assert False
    except ValueError:
        pass


def test_jtt_backend_if_pyvolve_available():
    try:
        import pyvolve  # noqa: F401
    except ImportError:
        return
    backend = build_r0_backend("jtt")
    log_R0 = backend.log_mutation_rates(["ACDEFG"], max_seq_len=6)
    assert log_R0.shape == (1, 6, 20)
    assert torch.allclose(
        torch.logsumexp(log_R0[0, 0], dim=-1), torch.tensor(0.0), atol=1e-5
    )
