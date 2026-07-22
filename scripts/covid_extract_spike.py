#!/usr/bin/env python3
"""
Extract the Spike (S) gene CDS from raw whole-genome SARS-CoV-2 FASTAs.

The raw files in data/covid/train/*_covid_seqs.fasta are complete ~29.9kb
genomes (unlike H3N2's already-single-segment HA data), so TreeSBM needs one
extra step before grouping: reference-coordinate-align each genome with
nextclade, then slice out the Spike CDS by its known reference coordinates
(Wuhan-Hu-1 / NC_045512.2, CDS 21563-25384, 1-based inclusive -- the
standard annotation used throughout SARS-CoV-2 genomics). nextclade puts
every input sequence into that same reference coordinate frame, so a fixed
column slice recovers each sequence's own Spike region, gaps stripped.

One-time setup (cluster, treesbm env):
    conda install -n treesbm -c bioconda -c conda-forge nextclade
    nextclade dataset get --name sars-cov-2 --output-dir data/covid/nextclade_dataset

Usage:
    python scripts/covid_extract_spike.py data/covid/train/africa_covid_seqs.fasta
    # -> data/covid/train/africa_spike.fasta (ungapped in-frame nt CDS, id = accession)
"""

import argparse
import os
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).parent.parent

# Wuhan-Hu-1 (NC_045512.2) Spike CDS, 1-based inclusive -> 0-based half-open slice.
SPIKE_START, SPIKE_END = 21562, 25384
MIN_LEN = 900  # < ~30% of the gene -> alignment failed/truncated for that record


def run_nextclade(src: Path, dataset_dir: Path, out_dir: Path, jobs: int) -> Path:
    aligned = out_dir / f"{src.stem}_aligned_genome.fasta"
    if aligned.exists() and aligned.stat().st_size > 0:
        return aligned
    cmd = ["nextclade", "run", f"--input-dataset={dataset_dir}",
           f"--output-fasta={aligned}", f"--jobs={jobs}", str(src)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"nextclade failed:\n{result.stderr[-2000:]}")
    return aligned


def extract_spike(aligned: Path, out_fasta: Path) -> int:
    records = []
    for rec in SeqIO.parse(aligned, "fasta"):
        region = str(rec.seq)[SPIKE_START:SPIKE_END].replace("-", "").upper()
        if len(region) >= MIN_LEN:
            records.append(SeqRecord(Seq(region), id=rec.id, description=""))
    SeqIO.write(records, out_fasta, "fasta")
    return len(records)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_fasta")
    ap.add_argument("--dataset-dir", default="data/covid/nextclade_dataset")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    args = ap.parse_args()

    src = ROOT / args.raw_fasta
    dataset_dir = ROOT / args.dataset_dir
    out_fasta = src.with_name(src.stem.replace("_covid_seqs", "") + "_spike.fasta")

    print(f"[nextclade] aligning {src.name} to reference coordinates ({args.jobs} jobs)...")
    aligned = run_nextclade(src, dataset_dir, src.parent, args.jobs)

    n = extract_spike(aligned, out_fasta)
    print(f"[extract] {n} Spike CDS sequences -> {out_fasta}")


if __name__ == "__main__":
    main()
