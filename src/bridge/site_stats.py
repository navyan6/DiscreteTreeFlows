"""
Per-column site statistics for hotspot selection / analysis.

Design (MSA-select → tree-apply):
  1. Discover hotspot columns from the TRAIN MSA (aligned columns): e.g. fraction
     of sequences whose AA differs from the column consensus/modal residue
     (``msa_mut_freq``), or Shannon entropy.
  2. Apply the resulting binary mask on the *tree bridge* via existing
     ``mut_hotspot_mask`` loss weighting / force — not a flat-MSA generative path.

Tree-edge / root→leaf mutation rates are secondary diagnostics (lit comparison,
sanity), not the default train selector.
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import torch

from src.bridge.losses import select_mut_hotspots

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
GAP_CHARS = frozenset("-X*")


def _aa_lookup_table() -> np.ndarray:
    tbl = np.full(256, 255, dtype=np.int16)
    for i, aa in enumerate(AA_VOCAB):
        tbl[ord(aa)] = i
    return tbl


def _iter_batches(dataset) -> Iterable[dict]:
    for i in range(len(dataset)):
        yield dataset[i]


def _valid_aa(a: str) -> bool:
    return bool(a) and a not in GAP_CHARS and a in AA_VOCAB


def _accumulate_seqs(counts: np.ndarray, sequences: Iterable[str], max_seq_len: int, tbl: np.ndarray) -> None:
    for seq in sequences:
        if not seq:
            continue
        arr = tbl[np.frombuffer(seq[:max_seq_len].encode("latin-1"), dtype=np.uint8)]
        valid = arr < 20
        if valid.any():
            pos = np.nonzero(valid)[0]
            np.add.at(counts, (pos, arr[valid]), 1)


def compute_msa_column_mut_freq_bundle(
    dataset,
    max_seq_len: int,
    *,
    sequences: list[str] | None = None,
) -> dict:
    """
    Primary hotspot score + CSV fields.

    msa_mut_freq[j] = (# seqs with AA ≠ modal) / (# valid AAs at j). Gaps/X/* excluded.
    """
    counts = np.zeros((max_seq_len, 20), dtype=np.int64)
    tbl = _aa_lookup_table()

    if sequences is not None:
        _accumulate_seqs(counts, sequences, max_seq_len, tbl)
    else:
        for batch in _iter_batches(dataset):
            _accumulate_seqs(counts, batch["seqs"].values(), max_seq_len, tbl)

    total = counts.sum(axis=1)
    modal_idx = counts.argmax(axis=1)
    modal_count = counts[np.arange(max_seq_len), modal_idx]
    non_cons = total - modal_count
    freq = np.zeros(max_seq_len, dtype=np.float64)
    nonempty = total > 0
    freq[nonempty] = non_cons[nonempty] / total[nonempty]
    consensus = np.array(
        [AA_VOCAB[i] if total[j] > 0 else "-" for j, i in enumerate(modal_idx)],
        dtype=object,
    )
    return {
        "msa_mut_freq": freq,
        "n_valid": total,
        "n_non_consensus": non_cons,
        "consensus_aa": consensus,
        "aa_counts": counts,
    }


def compute_msa_column_mut_freq(
    dataset,
    max_seq_len: int,
    *,
    sequences: list[str] | None = None,
) -> torch.Tensor:
    """Primary hotspot score: float32 [L] MSA mut-freq (frac ≠ consensus)."""
    return torch.tensor(
        compute_msa_column_mut_freq_bundle(dataset, max_seq_len, sequences=sequences)[
            "msa_mut_freq"
        ],
        dtype=torch.float32,
    )


def compute_tree_edge_mut_freq(dataset, max_seq_len: int) -> torch.Tensor:
    """Secondary: parent→child AA change rate over TRAIN trees."""
    bundle = compute_tree_mut_freq_bundle(dataset, max_seq_len)
    return torch.tensor(bundle["edge_mut_freq"], dtype=torch.float32)


def compute_tree_root_leaf_mut_freq(dataset, max_seq_len: int) -> torch.Tensor:
    """Secondary: fraction of leaves ≠ root at each column (TRAIN trees)."""
    bundle = compute_tree_mut_freq_bundle(dataset, max_seq_len)
    return torch.tensor(bundle["root_leaf_freq"], dtype=torch.float32)


def compute_tree_mut_freq_bundle(dataset, max_seq_len: int) -> dict:
    """Both secondary tree diagnostics + raw counts (vectorized per edge)."""
    edge_changes = np.zeros(max_seq_len, dtype=np.int64)
    edge_total = np.zeros(max_seq_len, dtype=np.int64)
    rl_diffs = np.zeros(max_seq_len, dtype=np.int64)
    rl_total = np.zeros(max_seq_len, dtype=np.int64)
    tbl = _aa_lookup_table()

    def _pair_stats(a: str, b: str, chg: np.ndarray, tot: np.ndarray) -> None:
        la = min(len(a), len(b), max_seq_len)
        if la == 0:
            return
        ia = tbl[np.frombuffer(a[:la].encode("latin-1"), dtype=np.uint8)]
        ib = tbl[np.frombuffer(b[:la].encode("latin-1"), dtype=np.uint8)]
        ok = (ia < 20) & (ib < 20)
        if not ok.any():
            return
        tot[:la][ok] += 1
        chg[:la][ok & (ia != ib)] += 1

    for batch in _iter_batches(dataset):
        seqs = batch["seqs"]
        for parent, child in batch["edges"]:
            _pair_stats(seqs[parent], seqs[child], edge_changes, edge_total)

        node_ids = batch["node_ids"]
        root_id = node_ids[batch["root_index"]]
        root_seq = seqs[root_id]
        has_children = {p for p, _ in batch["edges"]}
        leaves = [nid for nid in node_ids if nid not in has_children]
        for lid in leaves:
            if lid == root_id:
                continue
            _pair_stats(root_seq, seqs[lid], rl_diffs, rl_total)

    return {
        "edge_mut_freq": edge_changes.astype(np.float64) / np.maximum(edge_total, 1),
        "edge_changes": edge_changes,
        "edge_total": edge_total,
        "root_leaf_freq": rl_diffs.astype(np.float64) / np.maximum(rl_total, 1),
        "root_leaf_diffs": rl_diffs,
        "root_leaf_total": rl_total,
    }


def hotspot_mask_from_scores(
    scores: torch.Tensor | np.ndarray,
    topk: int | None = None,
    frac: float | None = None,
) -> torch.Tensor:
    """Top-k / top-frac binary mask on any per-column score."""
    return select_mut_hotspots(scores, topk=topk, frac=frac)


def normalize_score_01(scores: torch.Tensor) -> torch.Tensor:
    s = torch.as_tensor(scores, dtype=torch.float32)
    lo = float(s.min())
    hi = float(s.max())
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return torch.zeros_like(s)
    return ((s - lo) / (hi - lo)).clamp(0.0, 1.0)
