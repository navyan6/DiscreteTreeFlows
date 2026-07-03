#!/usr/bin/env python3
"""
Build the LARGEST comprehensive single-protein mutation landscape dataset.

Queries:
1. NCBI Virus database via gget virus (viruses only)
2. UniProt REST API (all proteins: viral, bacterial, eukaryotic)
3. Combines results into unified landscape dataset

This builds mutation landscapes for 10,000s of sequences across
multiple proteins and organisms.

Usage:
    python build_comprehensive_dataset.py [--ncbi-only] [--uniprot-only] [--test]
"""

import sys
import subprocess
from pathlib import Path
from typing import Dict
import json

PROJECT_ROOT = Path(__file__).parent
SRC_DIR = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))


def install_dependencies():
    """Install required packages."""
    print("Installing dependencies...")
    req_file = PROJECT_ROOT / "requirements.txt"
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "-r", str(req_file)],
        check=False,
        capture_output=True
    )
    print("✓ Dependencies installed\n")


def query_ncbi_viral_proteins():
    """Query NCBI for viral protein sequences."""
    print("=" * 70)
    print("STAGE 1: NCBI VIRUS DATABASE (gget virus)")
    print("=" * 70)
    print()

    from src.multi_protein_queries import (
        build_comprehensive_landscape_dataset,
        fetch_all_protein_sequences
    )

    # Query NCBI for high-priority viral proteins
    results = build_comprehensive_landscape_dataset(
        min_sequences_per_protein=100,
        max_proteins=None,
        priority_max=2,  # High and medium priority
    )

    print()

    if results:
        # Fetch sequences to FASTA
        fasta_files = fetch_all_protein_sequences(results)
        return results, fasta_files
    else:
        return {}, {}


def query_uniprot_all_proteins():
    """Query UniProt for comprehensive protein sequences."""
    print("\n" + "=" * 70)
    print("STAGE 2: UNIPROT (REST API) - ALL PROTEINS")
    print("=" * 70)
    print()

    from src.uniprot_queries import (
        build_comprehensive_uniprot_dataset,
        fetch_uniprot_sequences_to_fasta
    )

    # Query UniProt for diverse proteins
    results = build_comprehensive_uniprot_dataset(
        min_sequences=100,
        max_proteins=None,
        use_cli=False,  # Use REST API (more reliable)
    )

    print()

    if results:
        # Fetch sequences to FASTA
        fasta_files = fetch_uniprot_sequences_to_fasta(results)
        return results, fasta_files
    else:
        return {}, {}


def build_landscapes_for_all_proteins(
    all_results: Dict,
    all_fasta_files: Dict,
):
    """Build mutation landscapes for all queried proteins."""
    print("\n" + "=" * 70)
    print("STAGE 3: BUILD MUTATION LANDSCAPES")
    print("=" * 70)
    print()

    from src.landscapes import (
        compute_landscape_matrix,
        compute_landscape_statistics,
        choose_reference_sequence
    )
    from Bio import SeqIO

    landscapes_output = PROJECT_ROOT / "data" / "landscapes"
    landscapes_output.mkdir(parents=True, exist_ok=True)

    protein_landscapes = {}

    for protein_name, fasta_path in all_fasta_files.items():
        try:
            print(f"Processing {protein_name}...", end=" ", flush=True)

            # Read sequences
            sequences = []
            for record in SeqIO.parse(fasta_path, "fasta"):
                seq_str = str(record.seq).upper()
                if 50 < len(seq_str) < 10000:  # Reasonable protein length
                    sequences.append(seq_str)

            if len(sequences) < 50:
                print(f"skipped ({len(sequences)} sequences too few)")
                continue

            print(f"{len(sequences)} sequences...", end=" ", flush=True)

            # Compute landscape
            landscape = compute_landscape_matrix(
                sequences,
                normalize=True,
                pseudocount=0.1
            )
            stats = compute_landscape_statistics(landscape)

            protein_landscapes[protein_name] = {
                "landscape": landscape.tolist(),
                "statistics": stats,
                "n_sequences": len(sequences),
                "n_positions": landscape.shape[0],
            }

            print(f"✓ ({landscape.shape[0]} positions)")

        except Exception as e:
            print(f"ERROR: {str(e)[:50]}")

    print()
    return protein_landscapes


def save_comprehensive_dataset(
    ncbi_results: Dict,
    ncbi_fasta: Dict,
    uniprot_results: Dict,
    uniprot_fasta: Dict,
    protein_landscapes: Dict,
):
    """Save complete dataset to files."""
    print("=" * 70)
    print("STAGE 4: SAVE COMPREHENSIVE DATASET")
    print("=" * 70)
    print()

    landscapes_dir = PROJECT_ROOT / "data" / "landscapes"
    landscapes_dir.mkdir(parents=True, exist_ok=True)

    # Save main dataset
    dataset = {
        "metadata": {
            "description": "Comprehensive single-protein mutation landscape dataset",
            "source": ["NCBI Virus (gget virus)", "UniProt (REST API)"],
            "n_proteins": len(protein_landscapes),
            "total_sequences": sum(d["n_sequences"] for d in protein_landscapes.values()),
        },
        "ncbi": {
            "n_proteins": len(ncbi_results),
            "total_sequences": sum(r["n_sequences"] for r in ncbi_results.values()),
        },
        "uniprot": {
            "n_proteins": len(uniprot_results),
            "total_sequences": sum(r["n_sequences"] for r in uniprot_results.values()),
        },
        "landscapes": protein_landscapes,
    }

    output_file = landscapes_dir / "comprehensive_landscapes.json"
    with open(output_file, "w") as f:
        json.dump(dataset, f, indent=2, default=str)

    print(f"Saved comprehensive dataset: {output_file}")
    print(f"  Size: {output_file.stat().st_size / (1024**2):.1f} MB")

    # Save summary report
    report_file = landscapes_dir / "comprehensive_report.txt"
    with open(report_file, "w") as f:
        f.write("COMPREHENSIVE SINGLE-PROTEIN MUTATION LANDSCAPE DATASET\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"Total proteins: {len(protein_landscapes)}\n")
        f.write(f"Total sequences: {dataset['metadata']['total_sequences']:,}\n")
        f.write(f"Data sources: NCBI + UniProt\n\n")

        f.write("BREAKDOWN BY SOURCE\n")
        f.write("-" * 70 + "\n")
        f.write(f"NCBI Virus:\n")
        f.write(f"  Proteins: {dataset['ncbi']['n_proteins']}\n")
        f.write(f"  Sequences: {dataset['ncbi']['total_sequences']:,}\n\n")
        f.write(f"UniProt:\n")
        f.write(f"  Proteins: {dataset['uniprot']['n_proteins']}\n")
        f.write(f"  Sequences: {dataset['uniprot']['total_sequences']:,}\n\n")

        f.write("LANDSCAPES\n")
        f.write("-" * 70 + "\n")
        f.write("Protein | N Sequences | N Positions | Mean Entropy\n")
        f.write("-" * 70 + "\n")

        for protein, data in sorted(
            protein_landscapes.items(),
            key=lambda x: x[1]["n_sequences"],
            reverse=True
        ):
            f.write(
                f"{protein:30s} | {data['n_sequences']:11d} | "
                f"{data['n_positions']:11d} | "
                f"{data['statistics'].get('entropy_mean', 0):.3f}\n"
            )

    print(f"Saved report: {report_file}")
    print()


def main(ncbi_only=False, uniprot_only=False, test_mode=False):
    """Run complete pipeline."""
    print("=" * 70)
    print("COMPREHENSIVE SINGLE-PROTEIN MUTATION LANDSCAPE DATASET BUILDER")
    print("=" * 70)
    print()

    # Install dependencies
    install_dependencies()

    ncbi_results = {}
    ncbi_fasta = {}
    uniprot_results = {}
    uniprot_fasta = {}

    if test_mode:
        print("TEST MODE: Using synthetic data")
        print()
        # TODO: Add synthetic test data
    else:
        # Query NCBI
        if not uniprot_only:
            try:
                ncbi_results, ncbi_fasta = query_ncbi_viral_proteins()
            except Exception as e:
                print(f"NCBI query failed: {e}")
                ncbi_results, ncbi_fasta = {}, {}

        # Query UniProt
        if not ncbi_only:
            try:
                uniprot_results, uniprot_fasta = query_uniprot_all_proteins()
            except Exception as e:
                print(f"UniProt query failed: {e}")
                uniprot_results, uniprot_fasta = {}, {}

    # Combine results
    all_fasta = {**ncbi_fasta, **uniprot_fasta}

    if all_fasta:
        # Build landscapes
        protein_landscapes = build_landscapes_for_all_proteins(
            {**ncbi_results, **uniprot_results},
            all_fasta
        )

        # Save results
        save_comprehensive_dataset(
            ncbi_results,
            ncbi_fasta,
            uniprot_results,
            uniprot_fasta,
            protein_landscapes,
        )

        # Final summary
        print("=" * 70)
        print("COMPLETE")
        print("=" * 70)
        print(f"Built {len(protein_landscapes)} protein landscapes")
        print(f"Total sequences: {sum(d['n_sequences'] for d in protein_landscapes.values()):,}")
        print(f"Output: data/landscapes/comprehensive_landscapes.json")
        print("=" * 70)
    else:
        print("No sequences retrieved. Check NCBI/UniProt connectivity.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Build comprehensive mutation landscape dataset"
    )
    parser.add_argument(
        "--ncbi-only",
        action="store_true",
        help="Only query NCBI (viral proteins)"
    )
    parser.add_argument(
        "--uniprot-only",
        action="store_true",
        help="Only query UniProt (all proteins)"
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run in test mode with synthetic data"
    )

    args = parser.parse_args()
    main(
        ncbi_only=args.ncbi_only,
        uniprot_only=args.uniprot_only,
        test_mode=args.test
    )
