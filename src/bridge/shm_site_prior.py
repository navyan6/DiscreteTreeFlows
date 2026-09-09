"""
Antibody SHM-like site-rate prior on top of ESM Q0 (Recipe A).

Does **not** replace destination AA ratios (those stay ESM). Adjusts stay logits:
  - hot sites (CDR mask ∪ optional AID WRCH/DGYW): stay logit − boost
  - framework sites: stay logit + fwr_stay

AID motifs are detected on a reverse-translated NT string (most-common human
codons). Motifs: WRCH / DGYW (and RGYW/WRCY class). The middle C of WRCH and
the middle G of DGYW/RGYW mark the targeted AA column.

β fitness tilt should stay 0 for Abs when using this prior.
"""

from __future__ import annotations

from typing import Optional, Sequence, Union

import torch

from src.r0_backends import AA_TO_IDX

# Most-common human codon per AA (standard genetic code; stops unused).
_AA_TO_CODON = {
    "A": "GCC",
    "C": "TGC",
    "D": "GAC",
    "E": "GAG",
    "F": "TTC",
    "G": "GGC",
    "H": "CAC",
    "I": "ATC",
    "K": "AAG",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCC",
    "Q": "CAG",
    "R": "CGG",
    "S": "AGC",
    "T": "ACC",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAC",
}

# IUPAC: W=A/T, R=A/G, Y=C/T, D=A/G/T, H=A/C/T
_IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "W": "AT",
    "R": "AG",
    "Y": "CT",
    "D": "AGT",
    "H": "ACT",
    "N": "ACGT",
}


def _match_iupac(seq4: str, motif: str) -> bool:
    if len(seq4) != len(motif):
        return False
    for s, m in zip(seq4, motif):
        allowed = _IUPAC.get(m, m)
        if s not in allowed:
            return False
    return True


def reverse_translate_aa(aa: str) -> str:
    """AA string → NT via most-common codon; non-AA → NNN."""
    out = []
    for a in aa.upper():
        out.append(_AA_TO_CODON.get(a, "NNN"))
    return "".join(out)


def aid_hot_mask_from_aa(
    aa: str,
    max_seq_len: int,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.bool,
) -> torch.Tensor:
    """
    Boolean [L] mask: True where reverse-translated AID motifs hit that AA column.

    Motifs scanned on NT: WRCH (target = motif[2] = C → AA of that codon),
    DGYW / RGYW (target = motif[1] = G).
    """
    L = min(len(aa), max_seq_len)
    mask = torch.zeros(max_seq_len, dtype=dtype, device=device)
    if L <= 0:
        return mask
    nt = reverse_translate_aa(aa[:L])
    n = len(nt)
    # (motif, target_offset_within_motif)
    motifs = (("WRCH", 2), ("DGYW", 1), ("RGYW", 1), ("WRCY", 2))
    for i in range(max(0, n - 3)):
        window = nt[i : i + 4]
        if "N" in window:
            continue
        for motif, t_off in motifs:
            if not _match_iupac(window, motif):
                continue
            aa_idx = (i + t_off) // 3
            if 0 <= aa_idx < L:
                mask[aa_idx] = True
    return mask


def combine_hot_mask(
    cdr_mask: Optional[torch.Tensor],
    sequences: Optional[Sequence[str]],
    max_seq_len: int,
    use_aid: bool,
    device: torch.device,
) -> Optional[torch.Tensor]:
    """
    Return [L] bool hot mask, or [N, L] if AID is sequence-dependent and N>1.

    CDR mask is shared across the batch; AID may differ per sequence.
    """
    base: Optional[torch.Tensor] = None
    if cdr_mask is not None:
        base = cdr_mask.to(device=device, dtype=torch.bool).view(-1)
        if base.numel() < max_seq_len:
            pad = torch.zeros(max_seq_len - base.numel(), dtype=torch.bool, device=device)
            base = torch.cat([base, pad], dim=0)
        elif base.numel() > max_seq_len:
            base = base[:max_seq_len]

    if not use_aid or not sequences:
        return base

    aids = [
        aid_hot_mask_from_aa(seq, max_seq_len, device=device) for seq in sequences
    ]
    aid_stack = torch.stack(aids, dim=0)  # [N, L]
    if base is None:
        return aid_stack
    return aid_stack | base.view(1, -1)


def apply_shm_site_prior(
    log_R0: torch.Tensor,
    sequences: Sequence[str],
    *,
    cdr_mask: Optional[torch.Tensor] = None,
    boost: float = 2.0,
    fwr_stay: float = 0.5,
    use_aid: bool = True,
    enabled: bool = True,
) -> torch.Tensor:
    """
    Adjust stay logits on hot vs FWR sites.

    Args:
        log_R0: [N, L, 20] (or [L, 20] with N=1 sequences).
        sequences: length-N AA strings (parent / active leaf).
        cdr_mask: optional [L] CDR hotspot mask.
        boost: subtract from stay logit on hot sites (encourage mutation).
        fwr_stay: add to stay logit on non-hot sites (encourage conservation).
        use_aid: OR CDR with reverse-translated AID motif columns.
        enabled: if False, return log_R0 unchanged.
    """
    if not enabled or (boost == 0.0 and fwr_stay == 0.0 and not use_aid and cdr_mask is None):
        return log_R0
    if boost == 0.0 and fwr_stay == 0.0 and cdr_mask is None and not use_aid:
        return log_R0

    x = log_R0
    squeeze = False
    if x.ndim == 2:
        x = x.unsqueeze(0)
        squeeze = True
    if x.ndim != 3 or x.size(-1) != 20:
        raise ValueError(f"log_R0 must be [N,L,20], got {tuple(log_R0.shape)}")

    n, L, _ = x.shape
    if len(sequences) != n:
        raise ValueError(f"len(sequences)={len(sequences)} != batch N={n}")

    hot = combine_hot_mask(cdr_mask, sequences, L, use_aid, x.device)
    if hot is None:
        return log_R0
    if hot.ndim == 1:
        hot = hot.view(1, L).expand(n, L)
    else:
        hot = hot[:, :L]

    out = x.clone()
    for i, seq in enumerate(sequences):
        seq_u = (seq or "").upper()
        Li = min(len(seq_u), L)
        for j in range(Li):
            aa = seq_u[j]
            if aa not in AA_TO_IDX:
                continue
            stay = AA_TO_IDX[aa]
            if bool(hot[i, j]):
                if boost != 0.0:
                    out[i, j, stay] = out[i, j, stay] - float(boost)
            else:
                if fwr_stay != 0.0:
                    out[i, j, stay] = out[i, j, stay] + float(fwr_stay)
    return out.squeeze(0) if squeeze else out


def shm_prior_active(
    boost: float = 0.0,
    fwr_stay: float = 0.0,
    use_aid: bool = False,
    cdr_mask: Optional[torch.Tensor] = None,
) -> bool:
    return bool(
        (boost != 0.0)
        or (fwr_stay != 0.0)
        or use_aid
        or (cdr_mask is not None and float(boost) == 0.0 and float(fwr_stay) != 0.0)
        or (cdr_mask is not None and (boost != 0.0 or fwr_stay != 0.0))
    )
