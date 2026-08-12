"""Unit tests for COVID spike domain / region metric helpers."""

import json
from pathlib import Path

import torch

from src.spike_domains import (
    PRIMARY_DOMAIN_KEYS,
    REGION_RECALL_KEYS,
    SPIKE_DOMAINS_1BASED,
    all_domain_masks,
    domain_cols_0based,
    domain_mut_frac,
    region_any_mut,
    region_site_recall,
    spike_region_annotations,
    write_spike_region_annotations,
)


def test_rbd_band_matches_coord_validation():
    # COORD_VALIDATION: RBD_cols_0based [318, 540] ⇔ residues 319–541
    cols = list(domain_cols_0based(319, 541, 1280))
    assert cols[0] == 318
    assert cols[-1] == 540
    assert SPIKE_DOMAINS_1BASED["RBD"] == (319, 541)
    assert SPIKE_DOMAINS_1BASED["NTD"] == (13, 305)


def test_furin_s1_s2_boundary():
    assert SPIKE_DOMAINS_1BASED["furin"] == (680, 685)
    assert SPIKE_DOMAINS_1BASED["S1"][1] == 685
    assert SPIKE_DOMAINS_1BASED["S2"][0] == 686
    cols = list(domain_cols_0based(680, 685, 1273))
    assert cols == list(range(679, 685))  # 0-based 679..684


def test_domain_mut_frac_rbd_vs_ntd():
    L = 1280
    masks = all_domain_masks(L)
    # Mutations at N501 (col 500) and A222 (col 221)
    mut_cols = [500, 221]
    assert abs(domain_mut_frac(mut_cols, masks["RBD"]) - 0.5) < 1e-9
    assert abs(domain_mut_frac(mut_cols, masks["NTD"]) - 0.5) < 1e-9
    assert domain_mut_frac([], masks["RBD"]) != domain_mut_frac([], masks["RBD"])  # nan


def test_region_site_recall_and_any_mut():
    L = 600
    masks = all_domain_masks(L)
    # Build synthetic seqs: root all 'A'; GT mutates RBD col 400 and NTD col 100;
    # gen hits RBD only.
    root = "A" * L
    gt = list(root)
    gen = list(root)
    gt[100] = "V"   # NTD
    gt[400] = "Y"   # RBD
    gen[400] = "F"  # hits RBD site (wrong AA still counts)
    gt_s, gen_s = "".join(gt), "".join(gen)

    assert abs(region_site_recall(root, gt_s, gen_s, masks["RBD"]) - 1.0) < 1e-9
    assert abs(region_site_recall(root, gt_s, gen_s, masks["NTD"]) - 0.0) < 1e-9
    # No GT muts in furin → nan
    assert region_site_recall(root, gt_s, gen_s, masks["furin"]) != \
        region_site_recall(root, gt_s, gen_s, masks["furin"])

    assert region_any_mut(root, gen_s, masks["RBD"]) == 1.0
    assert region_any_mut(root, gen_s, masks["NTD"]) == 0.0
    assert region_any_mut(root, gt_s, masks["NTD"]) == 1.0


def test_primary_keys_have_masks():
    masks = all_domain_masks(1273)
    for k in PRIMARY_DOMAIN_KEYS:
        assert k in masks
        assert masks[k].dtype == torch.bool
        assert int(masks[k].sum()) > 0
    for k in REGION_RECALL_KEYS:
        assert k in masks


def test_write_spike_region_annotations(tmp_path: Path):
    out = write_spike_region_annotations(tmp_path / "region_annotations.json", L=1273)
    blob = json.loads(out.read_text())
    assert blob["pathogen"] == "covid_spike"
    assert "RBD" in blob["regions"]
    assert blob["regions"]["RBD"]["cols_0based"] == [318, 540]
    assert blob["regions"]["furin"]["lo_1based"] == 680
    assert "landmarks_1based" in blob["regions"]["RBD"]
    # smoke: spike_region_annotations shape
    ann = spike_region_annotations(1280)
    assert ann["regions"]["NTD"]["n_cols"] > 0
