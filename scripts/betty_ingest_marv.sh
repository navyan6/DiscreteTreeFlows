#!/bin/bash
# Append Marburg genomes, re-extract L, refresh outbreak splits.
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
cd ~/DiscreteTreeFlows

echo "=== Marburg NCBI ==="
$PYTHON scripts/download_ebolavirus_ncbi.py --species marv

echo "=== Re-extract L + window (pan-filo MSA) ==="
$PYTHON scripts/bdbv_extract_l.py
$PYTHON scripts/bdbv_define_l_window.py

echo "=== Outbreak splits ==="
$PYTHON scripts/prepare_filo_outbreak.py --out-base data/filo_l
$PYTHON scripts/audit_filo_splits.py --data data/filo_l

echo "Marburg ingest done: $(date)"
