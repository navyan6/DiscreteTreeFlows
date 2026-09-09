#!/usr/bin/env python3
"""Stratify aa|hit by polymerase-stall NT motifs (proxy reverse-translate).

TreeSBM is substitution-only AA. Homopolymer/G4/hairpin/repeat features live
on NT. We reverse-translate the root AA with most-common codon (same helper as
SHM AID prior) and tag each AA column. Then, for (root, GT leaf, nearest gen
leaf), split GT substitution sites into motif vs not.

This is a diagnostic proxy, not authentic genome NT.
"""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

from src.bridge.shm_site_prior import reverse_translate_aa

AA = set("ACDEFGHIKLMNPQRSTVWY")

# BLOSUM62 diagonal-ish pairs that count as "chemically close" misses.
BLOSUM_CLOSE = {
    ("D", "E"), ("D", "N"), ("E", "Q"), ("K", "R"), ("K", "Q"),
    ("I", "V"), ("I", "L"), ("I", "M"), ("L", "M"), ("L", "V"),
    ("F", "Y"), ("F", "W"), ("Y", "W"), ("S", "T"), ("A", "S"),
    ("N", "S"), ("H", "R"), ("H", "Q"),
}


def _close(a: str, b: str) -> bool:
    if a == b:
        return True
    p = (a, b) if a < b else (b, a)
    return p in BLOSUM_CLOSE


def parse_fasta(path: Path) -> dict[str, str]:
    recs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                recs[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        recs[name] = "".join(chunks)
    return recs


def homopolymer_aa_mask(nt: str, min_run: int = 4) -> list[bool]:
    L_aa = len(nt) // 3
    mask = [False] * L_aa
    i = 0
    n = len(nt)
    while i < n:
        j = i + 1
        while j < n and nt[j] == nt[i] and nt[i] in "ACGT":
            j += 1
        if j - i >= min_run:
            for k in range(i, j):
                mask[k // 3] = True
        i = j
    return mask


def g4_aa_mask(nt: str) -> list[bool]:
    L_aa = len(nt) // 3
    mask = [False] * L_aa
    for m in re.finditer(r"(G{3,}[ACGT]{1,7}){3,}G{3,}", nt):
        for k in range(m.start(), m.end()):
            mask[k // 3] = True
    return mask


def palindrome_aa_mask(nt: str, k: int = 8) -> list[bool]:
    """Inverted-repeat windows of length k (crude hairpin/palindrome proxy)."""
    L_aa = len(nt) // 3
    mask = [False] * L_aa
    comp = str.maketrans("ACGT", "TGCA")
    n = len(nt)
    for i in range(0, n - k + 1):
        w = nt[i : i + k]
        if "N" in w:
            continue
        rc = w.translate(comp)[::-1]
        # nearby reverse complement (gap 0..12) → stem-loop-ish
        for g in range(0, 13):
            j = i + k + g
            if j + k > n:
                break
            if nt[j : j + k] == rc:
                for t in range(i, j + k):
                    mask[t // 3] = True
                break
    return mask


def tandem_repeat_aa_mask(nt: str) -> list[bool]:
    L_aa = len(nt) // 3
    mask = [False] * L_aa
    for m in re.finditer(r"(.{2,6}?)\1{2,}", nt):
        for k in range(m.start(), m.end()):
            mask[k // 3] = True
    return mask


def motif_masks(aa: str) -> dict[str, list[bool]]:
    nt = reverse_translate_aa(aa)
    return {
        "homopolymer": homopolymer_aa_mask(nt),
        "g4": g4_aa_mask(nt),
        "palindrome": palindrome_aa_mask(nt),
        "tandem_repeat": tandem_repeat_aa_mask(nt),
    }


def hamming(a: str, b: str) -> int:
    n = min(len(a), len(b))
    return sum(x != y for x, y in zip(a[:n], b[:n])) + abs(len(a) - len(b))


def nearest(query: str, cands: list[str]) -> str:
    best, bd = cands[0], 10**9
    for c in cands:
        d = hamming(query, c)
        if d < bd:
            best, bd = c, d
    return best


def accumulate(root: str, gt: str, gen: str, masks: dict[str, list[bool]], bucket: dict):
    L = min(len(root), len(gt), len(gen))
    any_m = [any(masks[k][i] if i < len(masks[k]) else False for k in masks) for i in range(L)]
    for i in range(L):
        if root[i] not in AA or gt[i] not in AA or gen[i] not in AA:
            continue
        if root[i] == gt[i]:
            continue
        hit = gen[i] != root[i]
        correct = gen[i] == gt[i]
        close = _close(gen[i], gt[i]) if hit else False
        keys = ["all"]
        if any_m[i]:
            keys.append("any_stall")
        else:
            keys.append("no_stall")
        for name, m in masks.items():
            if i < len(m) and m[i]:
                keys.append(name)
        for k in keys:
            b = bucket[k]
            b["gt_mut"] += 1
            b["hit"] += int(hit)
            b["aa_ok"] += int(hit and correct)
            b["chem_ok"] += int(hit and close)


def summarize(bucket: dict) -> list[dict]:
    rows = []
    for k, b in sorted(bucket.items()):
        gt = b["gt_mut"] or 1
        hit = b["hit"]
        rows.append({
            "stratum": k,
            "gt_mut_sites": b["gt_mut"],
            "site_recall": b["hit"] / gt,
            "aa_hit": (b["aa_ok"] / hit) if hit else float("nan"),
            "chem_hit": (b["chem_ok"] / hit) if hit else float("nan"),
            "n_hits": hit,
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs", nargs="+", required=True)
    ap.add_argument("--gen", nargs="+", required=True)
    args = ap.parse_args()
    bucket = defaultdict(lambda: {"gt_mut": 0, "hit": 0, "aa_ok": 0, "chem_ok": 0})
    for op, gp in zip(args.obs, args.gen):
        obs = parse_fasta(Path(op))
        gen = parse_fasta(Path(gp))
        root = None
        for k, v in gen.items():
            if k.split("|")[0] == "root" or k == "root":
                root = v
                break
        if root is None:
            root = next(iter(obs.values()))
        masks = motif_masks(root)
        obs_leaves = [s for n, s in obs.items() if "root" not in n.lower()[:6]]
        gen_leaves = [s for n, s in gen.items() if not n.split("|")[0] in {"root"}]
        if not obs_leaves:
            obs_leaves = list(obs.values())[1:]
        if not gen_leaves:
            gen_leaves = list(gen.values())[1:]
        for gt in obs_leaves[:64]:
            g = nearest(gt, gen_leaves)
            accumulate(root, gt, g, masks, bucket)
    rows = summarize(bucket)
    print(f"{'stratum':<16} {'gt_mut':>8} {'n_hit':>8} {'site_rec':>9} {'aa|hit':>8} {'chem|hit':>9}")
    for r in rows:
        aa = r["aa_hit"]
        ch = r["chem_hit"]
        print(
            f"{r['stratum']:<16} {r['gt_mut_sites']:8d} {r['n_hits']:8d} "
            f"{r['site_recall']:9.3f} "
            f"{aa if aa==aa else float('nan'):8.3f} "
            f"{ch if ch==ch else float('nan'):9.3f}"
        )


if __name__ == "__main__":
    main()
