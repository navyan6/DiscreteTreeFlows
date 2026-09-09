#!/usr/bin/env python3
"""
Extract L-gene nucleotide CDS from filovirus whole-genome FASTAs.

Reads data/bdbv/manifest.json (or raw FASTA paths), aligns/slices L CDS by
reference coordinates per species, writes:
  data/bdbv/l_nt/{species}.fasta
  data/bdbv/l_nt/all.fasta

Headers: >ACC,YYYY-MM-DD,country,species,source
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scripts.bdbv_common import (
    DATA_ROOT,
    SPECIES,
    extract_l_nt_from_genbank,
    format_fasta_header,
    load_manifest,
    parse_fasta_header,
)

INSDC_RE = __import__("re").compile(r"^[A-Z]{1,2}\d+[A-Z]?\d*$")
PATHOPLEXUS_FASTA = DATA_ROOT / "raw" / "pathoplexus_2026" / "bdbv_genomes.fasta"


def _efetch_genbank(ids: list[str], retries: int = 3) -> list:
    import time

    from Bio import Entrez

    Entrez.email = "treesbm@example.com"
    records = []
    for i in range(0, len(ids), 200):
        chunk = ids[i : i + 200]
        for attempt in range(retries):
            try:
                handle = Entrez.efetch(
                    db="nucleotide", id=chunk, rettype="gb", retmode="text"
                )
                records.extend(SeqIO.parse(handle, "genbank"))
                handle.close()
                break
            except Exception as e:
                if attempt + 1 >= retries:
                    print(f"WARNING: GenBank efetch failed for chunk: {e}", flush=True)
                else:
                    time.sleep(1.5 * (attempt + 1))
        time.sleep(0.34)
    return records


def extract_from_manifest_genbank(
    manifest: list[dict],
    meta_by_acc: dict,
) -> dict[str, list[SeqRecord]]:
    by_species: dict[str, list[SeqRecord]] = {k: [] for k in SPECIES}
    for sp in SPECIES:
        accs = sorted(
            {
                m["accession"].split(".")[0]
                for m in manifest
                if m.get("species") == sp
                and m.get("accession")
                and INSDC_RE.match(m["accession"].split(".")[0])
            }
        )
        if not accs:
            continue
        print(f"[{sp}] GenBank refetch for {len(accs)} accessions", flush=True)
        for i in range(0, len(accs), 200):
            chunk = accs[i : i + 200]
            for rec in _efetch_genbank(chunk):
                acc = rec.id.split(".")[0]
                region = extract_l_nt_from_genbank(rec)
                if region is None:
                    region = _slice_l_from_genome(str(rec.seq).upper(), sp, try_all_species=False)
                if region is None:
                    continue
                m = meta_by_acc.get(acc, {})
                hdr = f"{m.get('date','unknown')},{m.get('country','unknown')},{m.get('species',sp)},{m.get('source','ncbi')}"
                by_species[sp].append(SeqRecord(Seq(region), id=acc, description=hdr))
            print(f"  [{sp}] {min(i+200, len(accs))}/{len(accs)} fetched, {len(by_species[sp])} L CDS", flush=True)
    return by_species


def _ref_genome_seq(species: str) -> str | None:
    ref_acc = SPECIES[species]["ref_accession"]
    ref_path = DATA_ROOT / "references" / f"{ref_acc}.fasta"
    if ref_path.is_file():
        return str(next(SeqIO.parse(ref_path, "fasta")).seq).upper()
    try:
        recs = _efetch_genbank([ref_acc.split(".")[0]])
        if recs:
            ref_path.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(recs[0], ref_path, "fasta")
            return str(recs[0].seq).upper()
    except Exception as e:
        print(f"WARNING: could not load reference {ref_acc}: {e}", flush=True)
    return None


def _slice_l_impute_from_ref(
    seq: str,
    ref_seq: str,
    species: str,
    min_resolved: float = 0.85,
) -> str | None:
    meta = SPECIES.get(species)
    if meta is None or len(seq) < meta["l_nt"][1] or len(ref_seq) < meta["l_nt"][1]:
        return None
    lo, hi = meta["l_nt"]
    region = str(seq[lo:hi]).upper()
    ref_region = str(ref_seq[lo:hi]).upper()
    if len(region) != len(ref_region) or len(region) % 3 != 0:
        return None
    resolved = sum(1 for c in region if c in "ACGT") / max(len(region), 1)
    if resolved < min_resolved:
        return None
    imputed = "".join(
        base if base in "ACGT" else ref_base
        for base, ref_base in zip(region, ref_region)
    )
    tr = str(Seq(imputed).translate())
    if "*" in tr[:-1]:
        return None
    return imputed


def _slice_l_from_genome(
    seq: str,
    species: str,
    try_all_species: bool = True,
    ref_seq: str | None = None,
    impute: bool = False,
) -> str | None:
    species_order = [species]
    if try_all_species:
        species_order += [k for k in SPECIES if k != species]

    for sp in species_order:
        meta = SPECIES.get(sp)
        if meta is None:
            continue
        lo, hi = meta["l_nt"]
        if len(seq) < hi:
            continue
        region = str(seq[lo:hi]).replace("-", "").upper()
        if len(region) < 5500 or len(region) % 3 != 0:
            continue
        tr = str(Seq(region).translate())
        if "*" in tr[:-1]:
            if impute and ref_seq is not None:
                imputed = _slice_l_impute_from_ref(seq, ref_seq, sp)
                if imputed is not None:
                    return imputed
            continue
        return region
    if impute and ref_seq is not None:
        return _slice_l_impute_from_ref(seq, ref_seq, species)
    return None


def _align_and_slice(ref_path: Path, query_seq: str, ref_lo: int, ref_hi: int) -> str | None:
    """Fallback: global align query to reference, map L coordinates."""
    ref_rec = next(SeqIO.parse(ref_path, "fasta"))
    ref_seq = str(ref_rec.seq).upper()
    tmp = DATA_ROOT / "_tmp_extract"
    tmp.mkdir(parents=True, exist_ok=True)
    q_f = tmp / "q.fasta"
    r_f = tmp / "r.fasta"
    SeqIO.write(SeqRecord(Seq(query_seq), id="q", description=""), q_f, "fasta")
    SeqIO.write(SeqRecord(Seq(ref_seq), id="r", description=""), r_f, "fasta")
    aln = tmp / "aln.fasta"
    try:
        subprocess.run(
            ["mafft", "--auto", "--quiet", str(r_f), str(q_f)],
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
        )
        # re-run with output file
        subprocess.run(
            ["mafft", "--auto", "--quiet", "-o", str(aln), str(r_f), str(q_f)],
            check=True,
            capture_output=True,
            timeout=120,
        )
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return _slice_l_from_genome(query_seq, "bdbv")

    aln_recs = list(SeqIO.parse(aln, "fasta"))
    if len(aln_recs) < 2:
        return None
    ref_aln = str(aln_recs[0].seq)
    q_aln = str(aln_recs[1].seq)
    # map ref ungapped positions to query
    ref_pos = 0
    q_chars: list[str] = []
    for a, b in zip(ref_aln, q_aln):
        if a == "-":
            continue
        ref_pos += 1
        if ref_lo < ref_pos <= ref_hi:
            if b != "-":
                q_chars.append(b)
    region = "".join(q_chars).upper()
    if len(region) < 6000 or len(region) % 3 != 0:
        return None
    return region


def extract_from_fasta(
    fasta_path: Path,
    species: str,
    records_meta: dict,
    align_fallback: bool = False,
    align_max: int = 50,
    impute_from_ref: bool = False,
) -> list[SeqRecord]:
    ref_path = DATA_ROOT / "references" / f"{SPECIES[species]['ref_accession']}.fasta"
    ref_seq = _ref_genome_seq(species) if impute_from_ref else None
    meta = SPECIES[species]
    out: list[SeqRecord] = []
    n_align = 0
    for i, rec in enumerate(SeqIO.parse(fasta_path, "fasta"), 1):
        acc = rec.id.split(",")[0].split()[0]
        seq = str(rec.seq).upper()
        region = _slice_l_from_genome(
            seq, species, ref_seq=ref_seq, impute=impute_from_ref
        )
        if region is None and align_fallback and ref_path.is_file() and n_align < align_max:
            region = _align_and_slice(ref_path, seq, meta["l_nt"][0], meta["l_nt"][1])
            n_align += 1
        if region is None:
            continue
        m = records_meta.get(acc, {})
        if isinstance(m, dict) and "date" in m:
            hdr = f"{m.get('date','unknown')},{m.get('country','unknown')},{m.get('species',species)},{m.get('source','ncbi')}"
        else:
            hdr = m.get("header_suffix") or f"unknown,unknown,{species},ncbi"
        out.append(SeqRecord(Seq(region), id=acc, description=hdr))
        if i % 100 == 0:
            print(f"  {fasta_path.name}: {i} scanned, {len(out)} L CDS", flush=True)
    if align_fallback and n_align:
        print(f"  {fasta_path.name}: MAFFT fallback used for {n_align} genomes", flush=True)
    return out


def _load_existing_l_nt(out_dir: Path) -> dict[str, list[SeqRecord]]:
    by_species: dict[str, list[SeqRecord]] = {k: [] for k in SPECIES}
    for sp in SPECIES:
        sp_path = out_dir / f"{sp}.fasta"
        if sp_path.is_file():
            by_species[sp] = list(SeqIO.parse(sp_path, "fasta"))
    return by_species


def _extract_pathoplexus_l(
    manifest: list[dict],
    meta_by_acc: dict,
    align_fallback: bool,
    align_max: int,
) -> list[SeqRecord]:
    out: list[SeqRecord] = []
    pp_manifest = [m for m in manifest if m.get("source") == "pathoplexus"]
    insdc_accs = sorted(
        {
            m["accession"].split(".")[0]
            for m in pp_manifest
            if INSDC_RE.match(m["accession"].split(".")[0])
        }
    )
    if insdc_accs:
        print(f"[pathoplexus] GenBank refetch for {len(insdc_accs)} INSDC accessions", flush=True)
        for rec in _efetch_genbank(insdc_accs):
            acc = rec.id.split(".")[0]
            region = extract_l_nt_from_genbank(rec) or _slice_l_from_genome(
                str(rec.seq).upper(), "bdbv", try_all_species=False
            )
            if region is None:
                continue
            m = meta_by_acc.get(acc, {})
            hdr = f"{m.get('date','unknown')},{m.get('country','unknown')},bdbv,pathoplexus"
            out.append(SeqRecord(Seq(region), id=acc, description=hdr))
    if PATHOPLEXUS_FASTA.is_file():
        out.extend(
            extract_from_fasta(
                PATHOPLEXUS_FASTA,
                "bdbv",
                meta_by_acc,
                align_fallback=align_fallback,
                align_max=align_max,
                impute_from_ref=True,
            )
        )
    seen: set[str] = set()
    uniq: list[SeqRecord] = []
    for rec in out:
        acc = rec.id.split()[0]
        if acc in seen:
            continue
        seen.add(acc)
        uniq.append(rec)
    return uniq


def _write_l_nt(
    by_species: dict[str, list[SeqRecord]],
    meta_by_acc: dict,
    out_dir: Path,
) -> int:
    all_recs: list[SeqRecord] = []
    for sp, recs in by_species.items():
        if not recs:
            continue
        seen: set[str] = set()
        uniq: list[SeqRecord] = []
        for r in recs:
            acc = r.id.split()[0].split(",")[0]
            if acc in seen:
                continue
            seen.add(acc)
            uniq.append(r)
        sp_path = out_dir / f"{sp}.fasta"
        with open(sp_path, "w") as f:
            for r in uniq:
                acc = r.id.split()[0].split(",")[0]
                m = meta_by_acc.get(acc, {})
                f.write(
                    format_fasta_header(
                        acc,
                        m.get("date", "unknown"),
                        m.get("country", "unknown"),
                        m.get("species", sp),
                        m.get("source", "ncbi"),
                    )
                    + f"\n{str(r.seq)}\n"
                )
        all_recs.extend(uniq)
        print(f"Wrote {len(uniq)} -> {sp_path}")
    all_path = out_dir / "all.fasta"
    with open(all_path, "w") as f:
        for r in all_recs:
            acc = r.id.split()[0].split(",")[0]
            m = meta_by_acc.get(acc, {})
            sp = m.get("species", "bdbv")
            f.write(
                format_fasta_header(
                    acc,
                    m.get("date", "unknown"),
                    m.get("country", "unknown"),
                    sp,
                    m.get("source", "ncbi"),
                )
                + f"\n{str(r.seq)}\n"
            )
    print(f"Total L CDS: {len(all_recs)} -> {all_path}")
    return len(all_recs)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default=str(DATA_ROOT / "manifest.json"))
    ap.add_argument(
        "--align-fallback",
        action="store_true",
        help="Use MAFFT align fallback for genomes that fail coordinate slice (slow)",
    )
    ap.add_argument("--align-max", type=int, default=50)
    ap.add_argument(
        "--fasta-only",
        action="store_true",
        help="Slice whole-genome FASTA only (fast but low yield; default uses GenBank L feature)",
    )
    ap.add_argument(
        "--pathoplexus-only",
        action="store_true",
        help="Merge Pathoplexus golden-test L into existing l_nt/ (no full manifest refetch)",
    )
    args = ap.parse_args()

    manifest = load_manifest()
    if not manifest:
        print("Empty manifest — run download_ebolavirus_ncbi.py first")
        sys.exit(1)

    meta_by_acc = {}
    for m in manifest:
        acc = m["accession"].split(".")[0]
        meta_by_acc[acc] = m

    out_dir = DATA_ROOT / "l_nt"
    out_dir.mkdir(parents=True, exist_ok=True)
    by_species: dict[str, list[SeqRecord]] = {k: [] for k in SPECIES}

    if args.pathoplexus_only:
        by_species = _load_existing_l_nt(out_dir)
        pp_accs = {
            m["accession"].split(".")[0]
            for m in manifest
            if m.get("source") == "pathoplexus"
        }
        for sp in SPECIES:
            by_species[sp] = [
                r for r in by_species.get(sp, []) if r.id.split()[0] not in pp_accs
            ]
        pp_recs = _extract_pathoplexus_l(
            manifest, meta_by_acc, args.align_fallback, args.align_max
        )
        by_species.setdefault("bdbv", []).extend(pp_recs)
        print(f"pathoplexus-only: merged {len(pp_recs)} golden-test L CDS")
        _write_l_nt(by_species, meta_by_acc, out_dir)
        return

    if args.fasta_only:
        fasta_paths = sorted(set(m.get("fasta", "") for m in manifest if m.get("fasta")))
        for rel in fasta_paths:
            fp = ROOT / rel
            if not fp.is_file():
                continue
            sp = fp.stem.split("_")[0]
            if sp not in SPECIES:
                for m in manifest:
                    if m.get("fasta") == rel:
                        sp = m["species"]
                        break
            recs = extract_from_fasta(
                fp, sp, meta_by_acc, align_fallback=args.align_fallback, align_max=args.align_max
            )
            by_species.setdefault(sp, []).extend(recs)
            print(f"{fp.name}: {len(recs)} L CDS")
    else:
        by_species = extract_from_manifest_genbank(manifest, meta_by_acc)
        for sp, recs in by_species.items():
            if recs:
                print(f"{sp}: {len(recs)} L CDS from GenBank features")
        if PATHOPLEXUS_FASTA.is_file():
            pp_recs = extract_from_fasta(
                PATHOPLEXUS_FASTA,
                "bdbv",
                meta_by_acc,
                align_fallback=args.align_fallback,
                align_max=args.align_max,
                impute_from_ref=True,
            )
            by_species.setdefault("bdbv", []).extend(pp_recs)
            print(f"pathoplexus: {len(pp_recs)} L CDS from {PATHOPLEXUS_FASTA.name}")

    _write_l_nt(by_species, meta_by_acc, out_dir)


if __name__ == "__main__":
    main()
