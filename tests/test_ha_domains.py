"""Unit tests for flu HA globular-head region helpers (H3 numbering)."""

import json
from pathlib import Path

import torch

from src.flu_lit_sites import (
    HA_DOMAINS_H3_1BASED,
    H3_GLOBULAR_HEAD_RANGE,
    H3_SIGNAL_LEN,
    PRIMARY_HA_DOMAIN_KEYS,
    all_ha_domain_masks,
    h3_to_col,
    ha_domain_cols_0based,
    ha_region_annotations,
    write_ha_region_annotations,
)
from src.spike_domains import domain_mut_frac, region_any_mut, region_site_recall


def test_ha_head_cols_match_coord_validation():
    # COORD_VALIDATION: globular_head_cols_0based [78, 267] for H3 #63–252
    lo, hi = H3_GLOBULAR_HEAD_RANGE
    assert (lo, hi) == (63, 252)
    cols = list(ha_domain_cols_0based(lo, hi, 566, H3_SIGNAL_LEN))
    assert cols[0] == 78
    assert cols[-1] == 267
    assert h3_to_col(63) == 78
    assert h3_to_col(252) == 267


def test_ha_head_region_metrics():
    L = 566
    masks = all_ha_domain_masks(L)
    assert set(masks) == set(PRIMARY_HA_DOMAIN_KEYS)
    m = masks["HA_head"]
    assert int(m.sum()) == 252 - 63 + 1

    root = "A" * L
    gt = list(root)
    gen = list(root)
    # F159 → col 174 (inside head); a stalk-ish col 400 (outside head)
    gt[174] = "S"
    gt[400] = "T"
    gen[174] = "Y"  # hit head site
    gt_s, gen_s = "".join(gt), "".join(gen)

    assert abs(region_site_recall(root, gt_s, gen_s, m) - 1.0) < 1e-9
    cols_gen = [174]
    assert abs(domain_mut_frac(cols_gen, m) - 1.0) < 1e-9
    assert region_any_mut(root, gen_s, m) == 1.0
    assert region_any_mut(root, "A" * L, m) == 0.0


def test_write_ha_region_annotations(tmp_path: Path):
    out = write_ha_region_annotations(tmp_path / "region_annotations.json", L=566)
    blob = json.loads(out.read_text())
    assert blob["pathogen"] == "flu_h3n2_ha"
    assert blob["regions"]["HA_head"]["cols_0based"] == [78, 267]
    assert blob["regions"]["HA_head"]["metric_prefix"] == "ha_head"
    ann = ha_region_annotations()
    assert HA_DOMAINS_H3_1BASED["HA_head"] == H3_GLOBULAR_HEAD_RANGE
    assert ann["indexing"]["signal_len"] == 16
    assert torch.is_tensor(all_ha_domain_masks(566)["HA_head"])
