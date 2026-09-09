"""Eval-time site-fire boost from polymerase-stall NT motifs.

Homopolymer / G4 / palindrome / tandem-repeat flags are computed on a
preferred-codon reverse-translate of the AA string (same helper as the SHM
AID prior). Flags scale *fire* (1 − p_stay), not destination mass.

Does not change checkpoints. Apply only at paint/eval.
"""

from __future__ import annotations

import re

import torch
import torch.nn.functional as F

from src.bridge.shm_site_prior import reverse_translate_aa
from src.r0_backends import AA_TO_IDX


def stall_site_mask(aa: str) -> list[bool]:
    nt = reverse_translate_aa(aa)
    L = len(aa)
    mask = [False] * L
    n = len(nt)

    i = 0
    while i < n:
        j = i + 1
        while j < n and nt[j] == nt[i] and nt[i] in "ACGT":
            j += 1
        if j - i >= 4:
            for k in range(i, j):
                mask[k // 3] = True
        i = j

    for m in re.finditer(r"(G{3,}[ACGT]{1,7}){3,}G{3,}", nt):
        for k in range(m.start(), m.end()):
            if k // 3 < L:
                mask[k // 3] = True
    for m in re.finditer(r"(.{2,6}?)\1{2,}", nt):
        for k in range(m.start(), m.end()):
            if k // 3 < L:
                mask[k // 3] = True

    comp = str.maketrans("ACGT", "TGCA")
    kmer = 8
    for i in range(0, n - kmer + 1):
        w = nt[i : i + kmer]
        if "N" in w:
            continue
        rc = w.translate(comp)[::-1]
        for g in range(0, 13):
            j = i + kmer + g
            if j + kmer > n:
                break
            if nt[j : j + kmer] == rc:
                for t in range(i, min(j + kmer, n)):
                    if t // 3 < L:
                        mask[t // 3] = True
                break
    return mask


def boost_fire_at_stalls(
    log_R: torch.Tensor,
    seq: str,
    scale: float,
) -> torch.Tensor:
    """log_R [L, 20] → same dest ratios, fire × (1+scale) on stall sites."""
    if scale <= 0:
        return log_R
    flags = stall_site_mask(seq)
    L = min(log_R.shape[0], len(seq), len(flags))
    probs = F.softmax(log_R[:L], dim=-1).clone()
    for pos in range(L):
        if not flags[pos]:
            continue
        a = AA_TO_IDX.get(seq[pos])
        if a is None:
            continue
        stay = float(probs[pos, a].item())
        fire = max(0.0, 1.0 - stay)
        fire2 = min(0.999, fire * (1.0 + float(scale)))
        stay2 = 1.0 - fire2
        off = probs[pos].clone()
        off[a] = 0.0
        s = float(off.sum().item())
        if s <= 0:
            continue
        probs[pos] = off / s * fire2
        probs[pos, a] = stay2
    out = log_R.clone()
    out[:L] = torch.log(probs.clamp_min(1e-12))
    return out
