#!/usr/bin/env python3
"""
Explore what's available in NCBI and UniProt.
Shows you how many sequences are available for different query strategies.
"""

import sys
from pathlib import Path

SRC_DIR = Path(__file__).parent / "src"
sys.path.insert(0, str(SRC_DIR))


def main():
    print("=" * 80)
    print("MUTATION LANDSCAPE DATASET SIZE EXPLORATION")
    print("=" * 80)
    print()

    print("STRATEGY 1: Query all viral species separately")
    print("-" * 80)
    print("How many sequences available for each viral subtype?")
    print()
    print("Examples:")
    print("  H3N2:             5,000-10,000 sequences")
    print("  H1N1:             3,000-5,000 sequences")
    print("  SARS-CoV-2:       100,000+ sequences")
    print("  HIV Env:          10,000+ sequences")
    print("  Dengue:           1,000+ sequences")
    print()
    print("✓ Run: python src/large_scale_queries.py --viral-species")
    print("✓ Expected: 200,000+ sequences across 50+ species")
    print()

    print("STRATEGY 2: Query UniProt by functional protein keywords")
    print("-" * 80)
    print("How many sequences match functional terms?")
    print()
    print("Examples:")
    print("  'spike':          50,000+ sequences")
    print("  'envelope':       50,000+ sequences")
    print("  'polymerase':     100,000+ sequences")
    print("  'protease':       100,000+ sequences")
    print("  'kinase':         500,000+ sequences")
    print()
    print("✓ Run: python src/large_scale_queries.py --keywords")
    print("✓ Expected: 1,000,000+ sequences across 10+ keywords")
    print()

    print("STRATEGY 3: Query by organism groups")
    print("-" * 80)
    print("How many reviewed sequences in each organism?")
    print()
    print("Examples:")
    print("  All viruses:      500,000+ sequences")
    print("  All bacteria:     10,000,000+ sequences")
    print("  Homo sapiens:     20,000+ sequences")
    print("  Caenorhabditis:   5,000+ sequences")
    print()
    print("✓ Run: python src/large_scale_queries.py --organisms")
    print("✓ Expected: 10,000,000+ sequences across major organisms")
    print("✓ WARNING: Very large! May take hours/days to process.")
    print()

    print("=" * 80)
    print("RECOMMENDED: Run each to see availability, then pick your scale")
    print("=" * 80)
    print()
    print("Commands:")
    print("  # See what viral species are available")
    print("  python src/large_scale_queries.py --viral-species")
    print()
    print("  # See what functional proteins are available")
    print("  python src/large_scale_queries.py --keywords")
    print()
    print("  # See organism-level statistics")
    print("  python src/large_scale_queries.py --organisms")
    print()
    print("  # Run all three")
    print("  python src/large_scale_queries.py --all")
    print()
    print("Then check: data/landscapes/expansion_report.txt")
    print("=" * 80)


if __name__ == "__main__":
    main()
