#!/usr/bin/env python3
"""Shared VOC threat-panel helpers (signatures, mut parsing, EVEscape lookup)."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

import torch
from Bio import SeqIO

AA_VOCAB = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_VOCAB)}
MUT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PANEL = ROOT / "benchmarks" / "voc_threat_panel.json"
DEFAULT_WT = ROOT / "data" / "covid" / "wt.txt"
DEFAULT_EVESCAPE = ROOT / "data" / "covid" / "evescape_spike_rbd.pt"


def load_panel(path: Path | str | None = None) -> dict:
    p = Path(path) if path else DEFAULT_PANEL
    return json.loads(p.read_text())


def voc_by_id(panel: dict, voc_id: str) -> dict:
    for v in panel["vocs"]:
        if v["id"] == voc_id:
            return v
    raise KeyError(f"VOC id not in panel: {voc_id}")


def parse_mut(mut: str) -> tuple[str, int, str]:
    """'K417T' -> (K, 417, T) with 1-based Spike position."""
    m = MUT_RE.match(mut.strip())
    if not m:
        raise ValueError(f"Bad mutation token: {mut!r}")
    return m.group(1), int(m.group(2)), m.group(3)


def load_wt_seq(path: Path | str | None = None) -> str:
    p = Path(path) if path else DEFAULT_WT
    rec = next(SeqIO.parse(p, "fasta"))
    return str(rec.seq).upper().replace("*", "")


_WT_CACHE: str | None = None
_MAP_CACHE: dict[str, tuple[int, ...]] = {}


def _default_wt() -> str:
    global _WT_CACHE
    if _WT_CACHE is None:
        _WT_CACHE = load_wt_seq() if DEFAULT_WT.exists() else ""
    return _WT_CACHE


def ref_pos_map(seq: str, ref: str | None = None) -> tuple[int, ...]:
    """Map each reference position (0-based) to its index in seq, or -1 if deleted.

    covid_extract_spike.py strips gaps after slicing Spike out of the
    reference-coordinate alignment, so a sequence carrying a Spike deletion is
    shorter than the reference and every residue downstream of the deletion sits
    at a lower index than its reference position. Indexing such a sequence
    directly by reference position silently reads the wrong residue -- that is
    what hid N501Y in Alpha and K417N in Beta (see SPIKE_FRAME_BUG.md).

    nextclade emits in reference coordinates, so insertions relative to the
    reference are dropped by the slice and length is exactly
    len(ref) - (number of deleted residues). Equal length therefore means no
    deletion and the identity map is correct, which is the common case and lets
    us skip alignment entirely.
    """
    ref = ref if ref is not None else _default_wt()
    if not ref or len(seq) == len(ref):
        return tuple(range(len(seq)))

    key = seq
    cached = _MAP_CACHE.get(key)
    if cached is not None:
        return cached

    from Bio import Align

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    mapping = [-1] * len(ref)
    try:
        aln = aligner.align(ref, seq)[0]
        for (rs, re_), (qs, _qe) in zip(*aln.aligned):
            for k in range(re_ - rs):
                mapping[rs + k] = qs + k
    except Exception:  # noqa: BLE001 - degenerate sequence; fall back to identity
        for i in range(min(len(ref), len(seq))):
            mapping[i] = i

    out = tuple(mapping)
    if len(_MAP_CACHE) < 50_000:
        _MAP_CACHE[key] = out
    return out


def seq_at(seq: str, pos1: int, ref: str | None = None) -> str | None:
    """1-based *reference* AA position -> AA on seq; None if deleted / OOB."""
    mapping = ref_pos_map(seq, ref)
    r = pos1 - 1
    if r < 0 or r >= len(mapping):
        return None
    i = mapping[r]
    if i < 0 or i >= len(seq):
        return None
    aa = seq[i]
    return aa if aa in AA_TO_IDX else None


def muts_present(seq: str, muts: Iterable[str], wt: str | None = None) -> list[str]:
    """Return signature muts whose mutant AA is present at the site on seq.

    If wt is given, also require seq differs from wt at that site (or matches
    the stated wt AA when checking consistency). We only require seq[pos]==mut.
    """
    hits = []
    for mut in muts:
        wt_aa, pos, mut_aa = parse_mut(mut)
        aa = seq_at(seq, pos)
        if aa == mut_aa:
            hits.append(mut)
        elif wt is not None and aa is not None:
            # optional consistency check — ignore mismatch with stated wt
            _ = wt_aa
    return hits


def muts_vs_root(seq: str, root: str, muts: Iterable[str]) -> list[str]:
    """Signature muts present on seq and different from root at that site."""
    hits = []
    for mut in muts:
        _, pos, mut_aa = parse_mut(mut)
        s = seq_at(seq, pos)
        r = seq_at(root, pos)
        if s == mut_aa and r is not None and s != r:
            hits.append(mut)
    return hits


def load_evescape(path: Path | str | None = None) -> tuple[torch.Tensor, dict]:
    p = Path(path) if path else DEFAULT_EVESCAPE
    blob = torch.load(p, map_location="cpu", weights_only=False)
    if isinstance(blob, dict):
        scores = blob["scores"]
        meta = {k: v for k, v in blob.items() if k != "scores"}
    else:
        scores = blob
        meta = {}
    return scores, meta


def evescape_of_mut(scores: torch.Tensor, mut: str) -> float | None:
    """Lookup EVEscape for a substitution; None if unscored / OOB."""
    _, pos, mut_aa = parse_mut(mut)
    i = pos - 1
    if i < 0 or i >= scores.shape[0]:
        return None
    j = AA_TO_IDX.get(mut_aa)
    if j is None:
        return None
    val = float(scores[i, j].item())
    return val if val != 0.0 else None  # 0 often = unscored in sparse RBD tensor


def annotate_muts_evescape(muts: Iterable[str], scores: torch.Tensor) -> dict[str, float | None]:
    return {m: evescape_of_mut(scores, m) for m in muts}


def load_fasta_seqs(path: Path | str) -> dict[str, str]:
    out = {}
    for rec in SeqIO.parse(path, "fasta"):
        out[rec.id] = str(rec.seq).upper().replace("-", "").replace("*", "")
    return out


def find_root_seq(seqs: dict[str, str], nwk_path: Path | None = None) -> tuple[str, str]:
    """Return (root_id, root_seq). Prefer NODE_0000001 / root|root / tree root name."""
    prefer = ["NODE_0000001", "root|root", "root"]
    for k in prefer:
        if k in seqs:
            return k, seqs[k]
    if nwk_path and Path(nwk_path).exists():
        try:
            from ete3 import Tree

            t = Tree(str(nwk_path), format=1)
            rname = t.name or (t.children[0].name if t.children else None)
            # ete rooted trees often have empty root name; use get_tree_root
            root = t.get_tree_root()
            # pick deepest internal with seq, else first leaf's ancestor with seq
            for node in [root] + list(root.traverse()):
                nm = node.name
                if nm and nm in seqs:
                    return nm, seqs[nm]
            if rname and rname in seqs:
                return rname, seqs[rname]
        except Exception:
            pass
    # fallback: shortest id starting with NODE_ or first key
    for k in sorted(seqs):
        if k.startswith("NODE_"):
            return k, seqs[k]
    k = next(iter(seqs))
    return k, seqs[k]


def leaf_seqs(seqs: dict[str, str]) -> dict[str, str]:
    """Heuristic: tips are accessions / |leaf suffixes, not NODE_/internal."""
    out = {}
    for k, v in seqs.items():
        if k.startswith("NODE_"):
            continue
        if k in ("root", "root|root") or k.endswith("|root"):
            continue
        if "|internal" in k:
            continue
        out[k] = v
    if not out:
        # generated trees use |leaf suffix
        out = {k: v for k, v in seqs.items() if k.endswith("|leaf") or "|leaf" in k}
    if not out:
        out = dict(seqs)
    return out


def score_group_against_voc(
    root_seq: str,
    leaves: dict[str, str],
    voc: dict,
) -> dict:
    sig = voc["signature"]
    bundle = voc["defining_bundle"]
    min_hits = int(voc.get("bundle_min_hits", len(bundle)))

    root_hits = muts_present(root_seq, sig)
    root_bundle = muts_present(root_seq, bundle)

    leaf_exact_counts = {m: 0 for m in sig}
    leaf_acquired_counts = {m: 0 for m in sig}
    n_bundle = 0
    n_bundle_acq = 0
    max_sig = 0
    max_acq = 0
    best_leaf = None

    for lid, seq in leaves.items():
        hits = muts_present(seq, sig)
        acq = muts_vs_root(seq, root_seq, sig)
        for m in hits:
            leaf_exact_counts[m] += 1
        for m in acq:
            leaf_acquired_counts[m] += 1
        bhits = muts_present(seq, bundle)
        bacq = muts_vs_root(seq, root_seq, bundle)
        if len(bhits) >= min_hits:
            n_bundle += 1
        if len(bacq) >= min_hits:
            n_bundle_acq += 1
        if len(hits) > max_sig:
            max_sig = len(hits)
            best_leaf = lid
        max_acq = max(max_acq, len(acq))

    n = max(len(leaves), 1)
    sites_hit = sum(1 for m, c in leaf_exact_counts.items() if c > 0)
    sites_acq = sum(1 for m, c in leaf_acquired_counts.items() if c > 0)

    return {
        "voc_id": voc["id"],
        "n_leaves": len(leaves),
        "n_signature": len(sig),
        "root_signature_hits": root_hits,
        "root_n_sig": len(root_hits),
        "root_bundle_hits": root_bundle,
        "root_has_full_bundle": len(root_bundle) >= min_hits,
        "prospective": len(root_bundle) < min_hits,
        "leaf_exact_counts": leaf_exact_counts,
        "leaf_acquired_counts": leaf_acquired_counts,
        "voc_exact_mut_recall": sites_hit / len(sig),
        "voc_acquired_mut_recall": sites_acq / len(sig),
        "voc_bundle_any": n_bundle > 0,
        "voc_bundle_frac": n_bundle / n,
        "voc_bundle_acquired_any": n_bundle_acq > 0,
        "voc_bundle_acquired_frac": n_bundle_acq / n,
        "max_sig_on_leaf": max_sig,
        "max_acquired_on_leaf": max_acq,
        "best_sig_leaf": best_leaf,
        "leaf_frac": {m: leaf_exact_counts[m] / n for m in sig},
        "acquired_frac": {m: leaf_acquired_counts[m] / n for m in sig},
    }


def group_paths(data_dir: Path, gid: int) -> dict[str, Path]:
    g3 = f"{gid:03d}"
    return {
        "anc": data_dir / f"group_{g3}_anc_aa.fasta",
        "nwk": data_dir / f"group_{g3}_rooted.nwk",
        "csv": data_dir / f"covidtest_group_{g3}.csv",
        "csv_train": data_dir / f"covidtrain_group_{g3}.csv",
        "csv_val": data_dir / f"covidval_group_{g3}.csv",
    }


def date_range_from_csv(csv_path: Path) -> tuple[str | None, str | None]:
    if not csv_path.exists():
        return None, None
    import csv as _csv

    dates = []
    with csv_path.open() as f:
        for row in _csv.DictReader(f):
            d = (row.get("date") or "").strip()
            if d and "X" not in d:
                dates.append(d)
    if not dates:
        return None, None
    return min(dates), max(dates)
