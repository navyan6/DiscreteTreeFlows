#!/bin/bash
#SBATCH --job-name=pv_inventory
#SBATCH --partition=genoa-std-mem
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --output=logs/panviral/inventory_%j.log

# Stage 1: count NCBI genomes for every eukaryotic virus species.
# Network-bound, not CPU-bound: ~53k taxa at the unauthenticated NCBI rate of
# 2.5 req/s is roughly 6 hours. Set NCBI_API_KEY to run at 9 req/s (~1.6 h).
# The counts cache makes this resumable, so a requeue costs nothing.

set -euo pipefail
REPO="${TREESBM_ROOT:-$HOME/DiscreteTreeFlows}"
cd "$REPO"
mkdir -p logs/panviral data/panviral

PY=${TREESBM_PY:-/vast/home/n/nnori/.conda/envs/treesbm/bin/python}

echo "host=$(hostname)  start=$(date -Is)  repo=$REPO"
$PY scripts/panviral/build_virus_inventory.py \
    --min-count "${MIN_COUNT:-150}" \
    --out data/panviral/virus_inventory.json
echo "done=$(date -Is)"
