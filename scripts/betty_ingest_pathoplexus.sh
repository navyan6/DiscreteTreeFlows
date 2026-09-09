#!/bin/bash
# Pathoplexus BDBV 2026 golden test ingest (eval-only; never train).
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
cd ~/DiscreteTreeFlows
mkdir -p logs data/bdbv/raw/pathoplexus_2026

echo "=== Pathoplexus Bdbv-2026 golden test (gz ok) ==="
$PYTHON scripts/bdbv_ingest_pathoplexus_export.py --search-dir data/bdbv

echo "Pathoplexus golden test ingest done: $(date)"
