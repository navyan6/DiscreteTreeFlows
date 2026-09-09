"""Sanity: Table 7 R0 / Table 8 lit-mask flags stay wired in CLIs."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_eval_evescape_has_r0_and_no_lit_mask():
    src = (ROOT / "scripts/eval_evescape_enrichment.py").read_text()
    assert "--r0-backend" in src
    assert "--no-lit-hotspot-mask" in src
    assert "r0_backend=r0_live" in src


def test_train_has_no_lit_hotspot_mask():
    src = (ROOT / "scripts/train.py").read_text()
    assert "--no-lit-hotspot-mask" in src
    assert "no_lit_hotspot_mask" in src


def test_generate_tree_accepts_r0_backend():
    src = (ROOT / "scripts/eval_single_tree.py").read_text()
    assert "r0_backend=None" in src
    assert "r0_backend.log_mutation_rates" in src


def test_table8_slurm_has_no_lit_mask_row():
    src = (ROOT / "scripts/slurm_table8_ablations.sh").read_text()
    assert "no_lit_mask" in src
    assert "--no-lit-hotspot-mask" in src
