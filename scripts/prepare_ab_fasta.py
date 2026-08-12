#!/usr/bin/env python3
"""
Parse / filter / shard data/Homo_sapiens.fasta for the antibody IgBLAST pipeline.

Input is a multi-source merge (~3.4M VH/VL-scale AA sequences; PLAbDab / SAbDab /
CoV-AbDab / IMGT / AbDb / PDB + BioProject PairedNGS dumps). Headers are messy
`|||`-joined annotations; this script extracts a flat TSV + clean FASTA shards
suitable for IgBLAST / Change-O clone clustering (Table 6 long pole).

Does NOT run IgBLAST, build clones, or write Newick trees.

Usage:
  # Full prep on Betty (CPU):
  python scripts/prepare_ab_fasta.py \\
    --input data/Homo_sapiens.fasta --out-dir data/ab_prep

  # Smoke / dry-run:
  python scripts/prepare_ab_fasta.py --input data/Homo_sapiens.fasta \\
    --out-dir data/ab_prep_smoke --max-seqs 5000 --n-shards 4

Outputs under --out-dir:
  shards/ab_shard_XXXX.fasta   clean >seq_id\\nSEQ (AA)
  metadata.tsv                 one row per kept sequence
  summary.json                 counts / length stats / source tallies
  clone_tag_counts.tsv         sparse _B-cells / clone-ish tags (if any)
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# Standard 20 AA + rare X / B / Z / J / U; gaps stripped later.
_AA_RE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWYXBZOJU]", re.I)
_VGENE_RE = re.compile(
    r"\b(IGH[VDJ][\w\*\-]+|IGK[VDJ][\w\*\-]+|IGL[VDJ][\w\*\-]+)\b", re.I
)
_CHAIN_HEAVY = re.compile(r"(heavy|IGH[VDJ]|_human_heavy|_heavy\b)", re.I)
_CHAIN_LIGHT = re.compile(r"(light|IGK[VDJ]|IGL[VDJ]|_human_light|_light\b)", re.I)
_BIOPROJECT_RE = re.compile(r"\b(PRJ[NE][A-Z]?\d+)\b", re.I)
_SRR_RE = re.compile(r"\b(SRR\d+)\b", re.I)
_CLONE_TAG_RE = re.compile(r"([A-Za-z0-9\-_]+_B-cells)", re.I)
_SOURCE_TAGS = (
    "PLAbDab",
    "SAbDab",
    "CoV-AbDab",
    "CoV-AbDab-PDB",
    "IMGT",
    "AbDb",
    "AbPDB",
    "PDB",
    "SACS",
    "PairedNGS",
)


def sanitize_aa(seq: str) -> str:
    s = seq.upper().replace(".", "").replace("-", "").replace("*", "")
    return _AA_RE.sub("", s)


def infer_chain(header: str, v_gene: str) -> str:
    """Best-effort heavy/light from trailing tags, then V gene, then keywords."""
    last = header.split("|||")[-1]
    if re.search(r"_heavy\b", last, re.I) or re.search(r"human_heavy", last, re.I):
        return "heavy"
    if re.search(r"_light\b", last, re.I) or re.search(r"human_light", last, re.I):
        return "light"
    if v_gene.startswith("IGH"):
        return "heavy"
    if v_gene.startswith(("IGK", "IGL")):
        return "light"
    primary = header.split("|||", 1)[0]
    if re.search(r"\bheavy\b", primary, re.I):
        return "heavy"
    if re.search(r"\blight\b", primary, re.I):
        return "light"
    if _CHAIN_HEAVY.search(header) and not _CHAIN_LIGHT.search(header):
        return "heavy"
    if _CHAIN_LIGHT.search(header) and not _CHAIN_HEAVY.search(header):
        return "light"
    return "unknown"


def parse_header(header: str) -> dict:
    """Extract best-effort fields from a multi-annotation FASTA header (no '>')."""
    h = header.strip()
    primary = h.split("|||", 1)[0].strip()
    seq_id = primary.split("|", 1)[0].split(";", 1)[0].strip() or "unknown"
    seq_id = re.sub(r"[^\w.\-]", "_", seq_id)[:120]

    v_gene = ""
    m = _VGENE_RE.search(h)
    if m:
        v_gene = m.group(1).upper().replace("_", "-")

    chain = infer_chain(h, v_gene)

    sources = []
    for tag in _SOURCE_TAGS:
        if tag in h:
            sources.append(tag)
    if _BIOPROJECT_RE.search(h):
        sources.append("BioProject")
    sources = sorted(set(sources)) or ["unknown"]

    bp = _BIOPROJECT_RE.search(h)
    srr = _SRR_RE.search(h)
    clone_m = _CLONE_TAG_RE.search(h)
    clone_tag = clone_m.group(1) if clone_m else ""

    return {
        "seq_id": seq_id,
        "v_gene": v_gene,
        "chain": chain,
        "source": "|".join(sources),
        "bioproject": bp.group(1) if bp else "",
        "srr": srr.group(1) if srr else "",
        "clone_tag": clone_tag,
        "header_short": primary[:200],
    }


def iter_fasta(path: Path, max_seqs: int | None = None):
    header = None
    seq_parts: list[str] = []
    n = 0
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                    n += 1
                    if max_seqs is not None and n >= max_seqs:
                        return
                header = line[1:].rstrip("\n")
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None and (max_seqs is None or n < max_seqs):
            yield header, "".join(seq_parts)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, default=ROOT / "data" / "Homo_sapiens.fasta")
    ap.add_argument("--out-dir", type=Path, default=ROOT / "data" / "ab_prep")
    ap.add_argument("--n-shards", type=int, default=64, help="Number of FASTA shards for IgBLAST")
    ap.add_argument("--min-len", type=int, default=80)
    ap.add_argument("--max-len", type=int, default=160, help="VH/VL-scale pad budget (Table 6)")
    ap.add_argument("--chains", default="heavy,light", help="Comma list: heavy,light,unknown")
    ap.add_argument("--dedup", action="store_true", default=True, help="Dedup by AA sequence (default on)")
    ap.add_argument("--no-dedup", action="store_false", dest="dedup")
    ap.add_argument("--prefer-heavy", action="store_true", help="If deduping, keep heavy over light for same AA")
    ap.add_argument("--max-seqs", type=int, default=None, help="Cap records read (smoke tests)")
    ap.add_argument("--progress-every", type=int, default=100_000)
    args = ap.parse_args()

    allow_chains = {c.strip().lower() for c in args.chains.split(",") if c.strip()}
    out = args.out_dir
    shard_dir = out / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)

    # Open shard handles lazily.
    shard_handles = []
    shard_counts = [0] * args.n_shards

    def shard_for(i: int):
        sid = i % args.n_shards
        if not shard_handles:
            for k in range(args.n_shards):
                p = shard_dir / f"ab_shard_{k:04d}.fasta"
                shard_handles.append(p.open("w"))
        return sid, shard_handles[sid]

    meta_path = out / "metadata.tsv"
    meta_f = meta_path.open("w", newline="")
    meta_w = csv.DictWriter(
        meta_f,
        fieldnames=[
            "seq_id",
            "uniq_id",
            "shard",
            "len",
            "v_gene",
            "chain",
            "source",
            "bioproject",
            "srr",
            "clone_tag",
            "header_short",
        ],
        delimiter="\t",
    )
    meta_w.writeheader()

    seen_seq: dict[str, str] = {}  # aa -> chain kept
    n_read = n_kept = n_filt_len = n_filt_chain = n_dedup = 0
    len_hist: Counter[int] = Counter()
    chain_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    v_family_counts: Counter[str] = Counter()
    clone_tag_counts: Counter[str] = Counter()

    if not args.input.is_file():
        raise SystemExit(f"Missing input FASTA: {args.input}")

    print(f"Reading {args.input} → {out}", flush=True)
    for header, raw_seq in iter_fasta(args.input, args.max_seqs):
        n_read += 1
        if args.progress_every and n_read % args.progress_every == 0:
            print(
                f"  … read={n_read:,} kept={n_kept:,} dedup_drop={n_dedup:,}",
                flush=True,
            )

        aa = sanitize_aa(raw_seq)
        L = len(aa)
        if L < args.min_len or L > args.max_len:
            n_filt_len += 1
            continue

        meta = parse_header(header)
        if meta["chain"] not in allow_chains:
            n_filt_chain += 1
            continue

        if args.dedup:
            prev = seen_seq.get(aa)
            if prev is not None:
                if args.prefer_heavy and prev != "heavy" and meta["chain"] == "heavy":
                    seen_seq[aa] = "heavy"
                    # Still skip writing a second copy; first record already written.
                n_dedup += 1
                continue
            seen_seq[aa] = meta["chain"]

        uniq_id = f"ab{n_kept:08d}_{meta['seq_id']}"[:80]
        sid, fh = shard_for(n_kept)
        fh.write(f">{uniq_id}\n")
        for i in range(0, L, 80):
            fh.write(aa[i : i + 80] + "\n")
        shard_counts[sid] += 1

        meta_w.writerow(
            {
                "seq_id": meta["seq_id"],
                "uniq_id": uniq_id,
                "shard": f"{sid:04d}",
                "len": L,
                "v_gene": meta["v_gene"],
                "chain": meta["chain"],
                "source": meta["source"],
                "bioproject": meta["bioproject"],
                "srr": meta["srr"],
                "clone_tag": meta["clone_tag"],
                "header_short": meta["header_short"],
            }
        )

        n_kept += 1
        len_hist[L] += 1
        chain_counts[meta["chain"]] += 1
        for s in meta["source"].split("|"):
            source_counts[s] += 1
        fam = re.match(r"(IGH[VDJ]\d+|IGK[VDJ]\d+|IGL[VDJ]\d+)", meta["v_gene"] or "")
        if fam:
            v_family_counts[fam.group(1)] += 1
        if meta["clone_tag"]:
            clone_tag_counts[meta["clone_tag"]] += 1

    for fh in shard_handles:
        fh.close()
    meta_f.close()

    # Length percentiles
    lengths = sorted(len_hist.elements())
    def pct(p: float) -> int | None:
        if not lengths:
            return None
        return lengths[min(len(lengths) - 1, int(p * (len(lengths) - 1)))]

    summary = {
        "input": str(args.input),
        "out_dir": str(out),
        "n_read": n_read,
        "n_kept": n_kept,
        "n_filt_len": n_filt_len,
        "n_filt_chain": n_filt_chain,
        "n_dedup_dropped": n_dedup,
        "min_len": args.min_len,
        "max_len": args.max_len,
        "n_shards": args.n_shards,
        "shard_counts": shard_counts,
        "chain_counts": dict(chain_counts),
        "source_counts": dict(source_counts.most_common()),
        "v_family_top": dict(v_family_counts.most_common(20)),
        "n_clone_tags": len(clone_tag_counts),
        "len_min": lengths[0] if lengths else None,
        "len_median": pct(0.5),
        "len_p90": pct(0.9),
        "len_max": lengths[-1] if lengths else None,
        "note": (
            "Shards ready for IgBLAST; clonal trees NOT built. "
            "See scripts/slurm_ab_igblast_clones.sh and ICLR_TABLE_FILL_PLAN §B."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    with (out / "clone_tag_counts.tsv").open("w", newline="") as cf:
        w = csv.writer(cf, delimiter="\t")
        w.writerow(["clone_tag", "n"])
        for tag, c in clone_tag_counts.most_common():
            w.writerow([tag, c])

    print(json.dumps({k: summary[k] for k in (
        "n_read", "n_kept", "n_filt_len", "n_filt_chain", "n_dedup_dropped",
        "len_min", "len_median", "len_p90", "len_max", "chain_counts",
    )}, indent=2))
    print(f"Wrote {meta_path} and {args.n_shards} shards under {shard_dir}", flush=True)


if __name__ == "__main__":
    main()
