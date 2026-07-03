#!/usr/bin/env python3
"""
Master phylogenetic pipeline: parse, QC, select, align, build trees, reconstruct.
Stage 1: Parse raw FASTA and quality-control the sequence pool.
"""
import argparse
import json
import os
from pathlib import Path
from collections import defaultdict, Counter
import re
from datetime import datetime

from Bio import SeqIO
import pandas as pd


def parse_gisaid_header(header):
    """Parse pipe-separated GISAID header into dict."""
    fields = header.split("|")
    keys = [
        "virus_name",
        "gisaid_isolate_id",
        "subtype_raw",
        "lineage",
        "clade",
        "collection_date",
        "submitter",
        "passage_history",
        "sample_provider_id",
        "submitting_lab_sample_id",
        "last_modified",
        "originating_lab",
        "submitting_lab",
        "segment",
        "segment_number",
        "segment_identifier",
        "dna_accession",
    ]

    parsed = {}
    for i, key in enumerate(keys):
        parsed[key] = fields[i] if i < len(fields) else ""

    # Extract subtype from A_/_H1N1 format
    subtype_raw = parsed["subtype_raw"]
    if "/" in subtype_raw:
        parsed["subtype"] = subtype_raw.split("/")[-1].strip("_")
    else:
        parsed["subtype"] = subtype_raw.strip("_")

    return parsed


def normalize_sequence(seq_str):
    """Uppercase, remove whitespace."""
    return seq_str.upper().strip()


def count_ambiguous(seq, valid_nucs="ACGT"):
    """Fraction of non-ACGT nucleotides."""
    n_invalid = sum(1 for c in seq if c not in valid_nucs)
    return n_invalid / len(seq) if seq else 0.0


def translate_to_protein(seq, frame=0):
    """Translate nucleotide sequence to protein."""
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }

    protein = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        aa = codon_table.get(codon, "X")
        protein.append(aa)

    return "".join(protein)


def has_internal_stops(protein):
    """Check for stop codons except at end."""
    return "*" in protein[:-1] if protein else False


def categorize_passage(passage_str):
    """Categorize passage history."""
    if not passage_str:
        return "unknown"
    s = passage_str.lower()
    if "original" in s:
        return "original"
    if "cell" in s and "egg" not in s:
        return "cell"
    if "egg" in s and "cell" not in s:
        return "egg"
    if "cell" in s and "egg" in s:
        return "mixed"
    return "unknown"


def stage_1_qc(
    input_fasta,
    output_dir,
    segment,
    subtype,
    lineage,
    length_bounds=None,
):
    """
    Stage 1: Parse and QC raw FASTA pool.
    """
    if length_bounds is None:
        length_bounds = {"HA": (1650, 1770), "NA": (1350, 1450)}

    os.makedirs(output_dir, exist_ok=True)

    parsed_records = []
    excluded_records = []
    qc_counts = defaultdict(int)

    # Read and parse
    for record in SeqIO.parse(input_fasta, "fasta"):
        header_dict = parse_gisaid_header(record.description)

        # Check segment and subtype
        if header_dict["segment"] != segment:
            excluded_records.append((
                record.id,
                "segment_mismatch",
                f"Expected {segment}, got {header_dict['segment']}",
            ))
            qc_counts["segment_mismatch"] += 1
            continue

        if header_dict["subtype"] != subtype:
            excluded_records.append((
                record.id,
                "subtype_mismatch",
                f"Expected {subtype}, got {header_dict['subtype']}",
            ))
            qc_counts["subtype_mismatch"] += 1
            continue

        if header_dict["lineage"] != lineage:
            excluded_records.append((
                record.id,
                "lineage_mismatch",
                f"Expected {lineage}, got {header_dict['lineage']}",
            ))
            qc_counts["lineage_mismatch"] += 1
            continue

        # Normalize sequence
        seq = normalize_sequence(str(record.seq))

        # Check length
        min_len, max_len = length_bounds.get(segment, (1200, 2000))
        if len(seq) < min_len or len(seq) > max_len:
            excluded_records.append((
                record.id,
                "length_out_of_bounds",
                f"Length {len(seq)}, bounds [{min_len}, {max_len}]",
            ))
            qc_counts["length_out_of_bounds"] += 1
            continue

        # Check ambiguous nucleotides
        amb_frac = count_ambiguous(seq)
        if amb_frac > 0.01:
            excluded_records.append((
                record.id,
                "high_ambiguity",
                f"Ambiguous fraction: {amb_frac:.3f}",
            ))
            qc_counts["high_ambiguity"] += 1
            continue

        # Check codons and internal stops
        protein = translate_to_protein(seq, frame=0)
        if has_internal_stops(protein):
            excluded_records.append((
                record.id,
                "internal_stop_codons",
                f"Protein: {protein}",
            ))
            qc_counts["internal_stop_codons"] += 1
            continue

        # Passage category
        passage_cat = categorize_passage(header_dict["passage_history"])

        row = {
            **header_dict,
            "record_id": record.id,
            "sequence_id": header_dict["segment_identifier"],
            "passage_category": passage_cat,
            "sequence_length": len(seq),
            "ambiguous_fraction": amb_frac,
            "protein_length": len(protein),
            "sequence": seq,
        }
        parsed_records.append(row)
        qc_counts["passed_qc"] += 1

    # Remove exact nucleotide duplicates within year/clade
    df_parsed = pd.DataFrame(parsed_records) if parsed_records else pd.DataFrame()

    duplicate_map = {}
    if len(df_parsed) > 0 and "record_id" in df_parsed.columns:
        df_parsed["year"] = pd.to_datetime(df_parsed["collection_date"]).dt.year
        dupes = df_parsed.groupby(["sequence", "year", "clade"]).size()
        dupes = dupes[dupes > 1]

        for (seq, year, clade), count in dupes.items():
            group = df_parsed[
                (df_parsed["sequence"] == seq) &
                (df_parsed["year"] == year) &
                (df_parsed["clade"] == clade)
            ]
            # Keep first, mark others as duplicates
            keep_id = group.iloc[0]["record_id"]
            for dup_id in group.iloc[1:]["record_id"]:
                duplicate_map[dup_id] = keep_id
                excluded_records.append((
                    dup_id,
                    "exact_duplicate",
                    f"Duplicate of {keep_id} in {year}/{clade}",
                ))
                qc_counts["exact_duplicate"] += 1

    # Remove duplicates from parsed
    df_filtered = df_parsed[~df_parsed["record_id"].isin(duplicate_map.keys())].copy()

    # Write outputs
    parsed_path = os.path.join(output_dir, "parsed_metadata.tsv")
    df_parsed[["record_id", "virus_name", "gisaid_isolate_id", "segment_identifier",
               "collection_date", "passage_category", "sequence_length", "ambiguous_fraction",
               "clade", "lineage", "subtype"]].to_csv(parsed_path, sep="\t", index=False)

    filtered_path = os.path.join(output_dir, "filtered_metadata.tsv")
    df_filtered[["record_id", "virus_name", "gisaid_isolate_id", "segment_identifier",
                 "collection_date", "passage_category", "sequence_length", "ambiguous_fraction",
                 "clade", "lineage", "subtype"]].to_csv(filtered_path, sep="\t", index=False)

    # Write filtered sequences
    fasta_path = os.path.join(output_dir, "filtered_nt.fasta")
    with open(fasta_path, "w") as f:
        for _, row in df_filtered.iterrows():
            f.write(f">{row['record_id']}\n{row['sequence']}\n")

    # Write excluded
    excluded_df = pd.DataFrame(excluded_records, columns=["record_id", "reason", "detail"])
    excluded_path = os.path.join(output_dir, "excluded_records.tsv")
    excluded_df.to_csv(excluded_path, sep="\t", index=False)

    # Write duplicate map
    dup_df = pd.DataFrame(list(duplicate_map.items()), columns=["duplicate_id", "kept_id"])
    dup_path = os.path.join(output_dir, "duplicate_map.tsv")
    dup_df.to_csv(dup_path, sep="\t", index=False)

    # Write QC summary
    summary = {
        "timestamp": datetime.now().isoformat(),
        "input_file": input_fasta,
        "segment": segment,
        "subtype": subtype,
        "lineage": lineage,
        "counts": dict(qc_counts),
        "total_input": len(df_parsed),
        "total_after_filtering": len(df_filtered),
        "total_excluded": len(excluded_records),
    }
    summary_path = os.path.join(output_dir, "qc_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*70}")
    print("STAGE 1: QUALITY CONTROL COMPLETE")
    print(f"{'='*70}")
    print(f"Input records:      {len(df_parsed)}")
    print(f"Passed QC:          {len(df_filtered)}")
    print(f"Excluded:           {len(excluded_records)}")
    print(f"\nQC breakdown:")
    for reason, count in sorted(qc_counts.items()):
        print(f"  {reason:30s} {count:6d}")
    print(f"\nOutputs in: {output_dir}")
    print(f"{'='*70}\n")

    return {
        "parsed_metadata": parsed_path,
        "filtered_metadata": filtered_path,
        "filtered_fasta": fasta_path,
        "excluded_records": excluded_path,
        "duplicate_map": dup_path,
        "qc_summary": summary_path,
    }


def stage_2_select(
    filtered_metadata_path,
    filtered_fasta_path,
    output_dir,
    target_size=1200,
    seed=2026,
):
    """
    Stage 2: Stratified selection of representatives.

    Stratify by year, then by clade (sqrt weighting), prefer original/cell passage,
    then random sample within each stratum.
    """
    import random
    random.seed(seed)

    os.makedirs(output_dir, exist_ok=True)

    # Load metadata
    df_meta = pd.read_csv(filtered_metadata_path, sep="\t")
    df_meta["year"] = pd.to_datetime(df_meta["collection_date"]).dt.year

    # Passage priority
    passage_priority = {"original": 0, "cell": 1, "unknown": 2, "egg": 3, "mixed": 4}
    df_meta["passage_rank"] = df_meta["passage_category"].map(
        lambda x: passage_priority.get(x, 5)
    )

    # Allocate targets per year
    year_counts = df_meta["year"].value_counts()
    years = sorted(year_counts.index)
    target_per_year = max(1, target_size // len(years))

    selected_ids = []

    for year in years:
        year_df = df_meta[df_meta["year"] == year].copy()
        year_target = min(target_per_year, len(year_df))

        # Allocate within year by clade (sqrt weighting)
        clade_counts = year_df["clade"].value_counts()
        clade_weights = {c: int(c_count ** 0.5) for c, c_count in clade_counts.items()}
        total_weight = sum(clade_weights.values())

        clade_targets = {}
        remaining = year_target
        for clade in sorted(clade_weights.keys()):
            target = max(1, int(year_target * clade_weights[clade] / total_weight))
            clade_targets[clade] = target
            remaining -= target

        if remaining > 0:
            clade_targets[sorted(clade_weights.keys())[0]] += remaining

        # Select per clade
        for clade, clade_target in clade_targets.items():
            clade_df = year_df[year_df["clade"] == clade].copy()
            if len(clade_df) == 0:
                continue

            # Sort by passage rank (prefer original)
            clade_df = clade_df.sort_values("passage_rank")

            if len(clade_df) <= clade_target:
                selected_ids.extend(clade_df["record_id"].tolist())
            else:
                # Random sample with passage preference
                selected_subset = random.sample(clade_df["record_id"].tolist(), clade_target)
                selected_ids.extend(selected_subset)

    # Ensure we have exactly target_size
    selected_ids = list(set(selected_ids))
    if len(selected_ids) < target_size:
        remaining_ids = [rid for rid in df_meta["record_id"] if rid not in selected_ids]
        shortfall = target_size - len(selected_ids)
        selected_ids.extend(random.sample(remaining_ids, min(shortfall, len(remaining_ids))))
    elif len(selected_ids) > target_size:
        selected_ids = random.sample(selected_ids, target_size)

    # Output selected sequences
    df_selected = df_meta[df_meta["record_id"].isin(selected_ids)].copy()

    selected_fasta = os.path.join(output_dir, "selected_1200_nt.fasta")
    with open(selected_fasta, "w") as f:
        for _, row in df_selected.iterrows():
            rec_id = row["record_id"]
            if rec_id in seq_dict:
                f.write(f">{rec_id}\n{seq_dict[rec_id]}\n")

    selected_meta = os.path.join(output_dir, "selected_1200_metadata.tsv")
    df_selected[["record_id", "virus_name", "collection_date", "passage_category",
                 "clade", "lineage", "year"]].to_csv(selected_meta, sep="\t", index=False)

    selected_accessions = os.path.join(output_dir, "selected_accessions.txt")
    with open(selected_accessions, "w") as f:
        for rec_id in sorted(selected_ids):
            f.write(f"{rec_id}\n")

    # Selection summary
    summary = {
        "target_size": target_size,
        "actual_size": len(selected_ids),
        "years_represented": len(years),
        "clades_represented": df_selected["clade"].nunique(),
        "passage_counts_before": dict(df_meta["passage_category"].value_counts()),
        "passage_counts_after": dict(df_selected["passage_category"].value_counts()),
    }

    summary_path = os.path.join(output_dir, "selection_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*70}")
    print("STAGE 2: REPRESENTATIVE SELECTION COMPLETE")
    print(f"{'='*70}")
    print(f"Target size:             {target_size}")
    print(f"Actual selected:         {len(selected_ids)}")
    print(f"Years represented:       {len(years)} ({min(years)}–{max(years)})")
    print(f"Clades represented:      {df_selected['clade'].nunique()}")
    print(f"\nPassage distribution after selection:")
    for passage, count in df_selected["passage_category"].value_counts().items():
        print(f"  {passage:15s} {count:6d}")
    print(f"\nOutputs in: {output_dir}")
    print(f"{'='*70}\n")

    return {
        "selected_fasta": selected_fasta,
        "selected_metadata": selected_meta,
        "selected_accessions": selected_accessions,
        "selection_summary": summary_path,
    }


def stage_3_align(
    filtered_fasta_path,
    output_dir,
):
    """
    Stage 3: MAFFT multiple sequence alignment.
    """
    import subprocess

    os.makedirs(output_dir, exist_ok=True)

    aligned_fasta = os.path.join(output_dir, "aligned_nt.fasta")

    # Run MAFFT
    cmd = f"mafft --auto {filtered_fasta_path} > {aligned_fasta}"
    print(f"\nRunning: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"MAFFT error: {result.stderr}")
        return None

    # Validate alignment
    aligned_records = list(SeqIO.parse(aligned_fasta, "fasta"))
    if not aligned_records:
        print("Error: No aligned sequences")
        return None

    align_len = len(aligned_records[0].seq)
    invalid = sum(1 for rec in aligned_records if len(rec.seq) != align_len)

    if invalid > 0:
        print(f"Warning: {invalid} sequences have different length")

    # Translate to protein
    aa_fasta = os.path.join(output_dir, "aligned_aa.fasta")
    with open(aa_fasta, "w") as f:
        for rec in aligned_records:
            protein = translate_to_protein(str(rec.seq), frame=0)
            f.write(f">{rec.id}\n{protein}\n")

    print(f"\n{'='*70}")
    print("STAGE 3: ALIGNMENT COMPLETE")
    print(f"{'='*70}")
    print(f"Sequences aligned:  {len(aligned_records)}")
    print(f"Alignment length:   {align_len} nt")
    print(f"Outputs in: {output_dir}")
    print(f"{'='*70}\n")

    return {
        "aligned_fasta": aligned_fasta,
        "aligned_aa": aa_fasta,
    }


def stage_4_iqtree(
    aligned_fasta_path,
    output_dir,
    seed=2026,
):
    """
    Stage 4: IQ-TREE maximum-likelihood tree.
    """
    import subprocess

    os.makedirs(output_dir, exist_ok=True)

    # Run IQ-TREE
    cmd = (
        f"iqtree3 "
        f"-s {aligned_fasta_path} "
        f"-m GTR+F+G4 "
        f"-B 1000 "
        f"-alrt 1000 "
        f"--seed {seed} "
        f"--prefix {os.path.join(output_dir, 'master')} "
        f"-quiet"
    )

    print(f"\nRunning IQ-TREE...")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"IQ-TREE error: {result.stderr}")
        return None

    treefile = os.path.join(output_dir, "master.treefile")
    if not os.path.exists(treefile):
        print("Error: IQ-TREE did not produce treefile")
        return None

    print(f"\n{'='*70}")
    print("STAGE 4: IQ-TREE COMPLETE")
    print(f"{'='*70}")
    print(f"Tree file:     {treefile}")
    print(f"Outputs in: {output_dir}")
    print(f"{'='*70}\n")

    return {
        "treefile": treefile,
        "contree": os.path.join(output_dir, "master.contree"),
        "iqtree_log": os.path.join(output_dir, "master.iqtree"),
    }


def stage_5_treetime(
    treefile_path,
    aligned_fasta_path,
    filtered_metadata_path,
    output_dir,
):
    """
    Stage 5: TreeTime temporal refinement with Augur.
    """
    import subprocess

    os.makedirs(output_dir, exist_ok=True)

    # Load and prepare metadata for Augur
    df_meta = pd.read_csv(filtered_metadata_path, sep="\t")

    # Minimal metadata: strain + date only (no extra columns)
    df_meta_aug = df_meta[["record_id", "collection_date"]].copy()
    df_meta_aug.columns = ["strain", "date"]

    # Save for Augur (only 2 columns)
    treetime_meta = os.path.join(output_dir, "metadata.tsv")
    df_meta_aug.to_csv(treetime_meta, sep="\t", index=False)

    # Run Augur refine
    refined_tree = os.path.join(output_dir, "refined_master.nwk")
    refined_nodes = os.path.join(output_dir, "refined_master_nodes.json")

    cmd = (
        f"augur refine "
        f"--tree {treefile_path} "
        f"--alignment {aligned_fasta_path} "
        f"--metadata {treetime_meta} "
        f"--timetree "
        f"--coalescent opt "
        f"--date-confidence "
        f"--output-tree {refined_tree} "
        f"--output-node-data {refined_nodes} "
        f"2>&1"
    )

    print(f"\nRunning Augur refine (TreeTime)...")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Augur refine error: {result.stderr}")
        return None

    if not os.path.exists(refined_tree):
        print("Error: TreeTime did not produce refined tree")
        return None

    print(f"\n{'='*70}")
    print("STAGE 5: TREETIME REFINEMENT COMPLETE")
    print(f"{'='*70}")
    print(f"Refined tree: {refined_tree}")
    print(f"Outputs in: {output_dir}")
    print(f"{'='*70}\n")

    return {
        "refined_tree": refined_tree,
        "refined_nodes": refined_nodes,
        "treetime_metadata": treetime_meta,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phylogenetic pipeline: parse, QC, select")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--segment", required=True, help="Segment (HA, NA, etc)")
    parser.add_argument("--subtype", required=True, help="Subtype (H1N1, H3N2, etc)")
    parser.add_argument("--lineage", required=True, help="Lineage (pdm09, etc)")
    parser.add_argument("--target-size", type=int, default=1200, help="Target # of tips")
    parser.add_argument("--seed", type=int, default=2026, help="Random seed")
    parser.add_argument("--stage", default="all", help="Stage to run (1, 2, or all)")

    args = parser.parse_args()

    base_output_dir = args.output_dir

    if args.stage in ["1", "all"]:
        stage1_dir = os.path.join(base_output_dir, "stage_1_qc")
        results_1 = stage_1_qc(
            input_fasta=args.input_fasta,
            output_dir=stage1_dir,
            segment=args.segment,
            subtype=args.subtype,
            lineage=args.lineage,
        )

    if args.stage in ["2", "all"]:
        stage1_dir = os.path.join(base_output_dir, "stage_1_qc")
        stage2_dir = os.path.join(base_output_dir, "stage_2_select")
        results_2 = stage_2_select(
            filtered_metadata_path=os.path.join(stage1_dir, "filtered_metadata.tsv"),
            filtered_fasta_path=os.path.join(stage1_dir, "filtered_nt.fasta"),
            output_dir=stage2_dir,
            target_size=args.target_size,
            seed=args.seed,
        )

    if args.stage in ["3", "all"]:
        stage1_dir = os.path.join(base_output_dir, "stage_1_qc")
        stage3_dir = os.path.join(base_output_dir, "stage_3_align")
        results_3 = stage_3_align(
            filtered_fasta_path=os.path.join(stage1_dir, "filtered_nt.fasta"),
            output_dir=stage3_dir,
        )

    if args.stage in ["4", "all"]:
        stage3_dir = os.path.join(base_output_dir, "stage_3_align")
        stage4_dir = os.path.join(base_output_dir, "stage_4_iqtree")
        results_4 = stage_4_iqtree(
            aligned_fasta_path=os.path.join(stage3_dir, "aligned_nt.fasta"),
            output_dir=stage4_dir,
            seed=args.seed,
        )

    if args.stage in ["5", "all"]:
        stage1_dir = os.path.join(base_output_dir, "stage_1_qc")
        stage3_dir = os.path.join(base_output_dir, "stage_3_align")
        stage4_dir = os.path.join(base_output_dir, "stage_4_iqtree")
        stage5_dir = os.path.join(base_output_dir, "stage_5_treetime")
        results_5 = stage_5_treetime(
            treefile_path=os.path.join(stage4_dir, "master.treefile"),
            aligned_fasta_path=os.path.join(stage3_dir, "aligned_nt.fasta"),
            filtered_metadata_path=os.path.join(stage1_dir, "filtered_metadata.tsv"),
            output_dir=stage5_dir,
        )
