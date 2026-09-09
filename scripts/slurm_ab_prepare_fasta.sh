#!/bin/bash
#SBATCH --job-name=ab_prep_fasta
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
# omit --mem: PARCC filter requires exact 5632MB/CPU on genoa-std-mem
#SBATCH --time=04:00:00
#SBATCH --output=logs/ab_prep_fasta_%j.log
#SBATCH --error=logs/ab_prep_fasta_%j.log
#
# Parse / filter / shard data/Homo_sapiens.fasta for IgBLAST (Table 6 long pole).
# CPU-only. Does NOT run IgBLAST or build trees.
#
# Usage (Betty login shell):
#   cd ~/DiscreteTreeFlows && mkdir -p logs
#   # Ensure data/Homo_sapiens.fasta is present (~1.5G), then:
#   sbatch scripts/slurm_ab_prepare_fasta.sh
#   # Overrides:
#   INPUT=data/Homo_sapiens.fasta OUT=data/ab_prep N_SHARDS=64 \
#     sbatch scripts/slurm_ab_prepare_fasta.sh
#   # Smoke:
#   MAX_SEQS=5000 N_SHARDS=4 OUT=data/ab_prep_smoke \
#     sbatch scripts/slurm_ab_prepare_fasta.sh

set -euo pipefail

INPUT="${INPUT:-data/Homo_sapiens.fasta}"
OUT="${OUT:-data/ab_prep}"
N_SHARDS="${N_SHARDS:-64}"
MIN_LEN="${MIN_LEN:-80}"
MAX_LEN="${MAX_LEN:-160}"
MAX_SEQS="${MAX_SEQS:-}"
CHAINS="${CHAINS:-heavy,light}"

mkdir -p logs "$OUT"

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows

echo "Start: $(date)"
echo "INPUT=$INPUT OUT=$OUT N_SHARDS=$N_SHARDS MIN_LEN=$MIN_LEN MAX_LEN=$MAX_LEN"

EXTRA=()
if [[ -n "$MAX_SEQS" ]]; then
  EXTRA+=(--max-seqs "$MAX_SEQS")
fi

$PYTHON scripts/prepare_ab_fasta.py \
  --input "$INPUT" \
  --out-dir "$OUT" \
  --n-shards "$N_SHARDS" \
  --min-len "$MIN_LEN" \
  --max-len "$MAX_LEN" \
  --chains "$CHAINS" \
  --prefer-heavy \
  "${EXTRA[@]}"

echo "Done: $(date)"
echo "Next: sbatch scripts/slurm_ab_igblast_clones.sh  (after IgBLAST+IMGT refs exist)"
