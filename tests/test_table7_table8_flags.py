"""Sanity: Table 7 R0 / Table 8 lit-mask flags stay wired in CLIs."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_eval_evescape_has_r0_and_no_lit_mask():
    src = (ROOT / "scripts/eval_evescape_enrichment.py").read_text()
    assert "--r0-backend" in src
    assert "--no-lit-hotspot-mask" in src
    assert "--ablate-site-entropy" in src
    assert "r0_backend=r0_live" in src
    assert "plm_nll" in src


def test_train_has_no_lit_hotspot_mask():
    src = (ROOT / "scripts/train.py").read_text()
    assert "--no-lit-hotspot-mask" in src
    assert "no_lit_hotspot_mask" in src


def test_generate_tree_accepts_r0_backend():
    src = (ROOT / "scripts/eval_single_tree.py").read_text()
    assert "r0_backend=None" in src
    assert "r0_backend.log_mutation_rates" in src
    assert "ablate_site_entropy" in src


def test_table8_slurm_has_no_lit_mask_row():
    src = (ROOT / "scripts/slurm_table8_ablations.sh").read_text()
    assert "no_lit_mask" in src
    assert "--no-lit-hotspot-mask" in src
    assert "no_entropy" in src
    assert "--ablate-site-entropy" in src
    assert 'MODE=baselines' in src or "baselines" in src


def test_progen2_not_stubbed():
    from src.r0_backends import STUB_BACKENDS, BACKEND_PROGEN2, BACKEND_EVO2
    assert BACKEND_PROGEN2 not in STUB_BACKENDS
    assert BACKEND_EVO2 in STUB_BACKENDS
    src = (ROOT / "src/r0_backends.py").read_text()
    assert "class ProGen2R0Backend" in src
