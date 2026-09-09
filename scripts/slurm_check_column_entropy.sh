#!/bin/bash
#SBATCH --job-name=col_entropy
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/col_entropy_%j.log
#SBATCH --error=logs/col_entropy_%j.log
#
# Column-entropy mut vs conserved sanity check (CPU-only; genoa-std-mem).
# No GPU / mig QOS needed. Override partition: sbatch -p <part> ...
#
# Usage:
#   sbatch scripts/slurm_check_column_entropy.sh <data_dir> [max_seq_len] [out_dir]
#   sbatch scripts/slurm_check_column_entropy.sh data/train 566 results/entropy_h3n2
#   sbatch scripts/slurm_check_column_entropy.sh data/covid/train 1280 results/entropy_covid

set -e

DATA_DIR="${1:?usage: sbatch scripts/slurm_check_column_entropy.sh <data_dir> [max_seq_len] [out_dir]}"
MAX_SEQ_LEN="${2:-566}"
OUT_DIR="${3:-results/column_entropy_sanity}"

mkdir -p logs "$OUT_DIR"

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows

echo "Start: $(date)"
echo "Data: $DATA_DIR  max_seq_len=$MAX_SEQ_LEN  out=$OUT_DIR"

$PYTHON scripts/check_column_entropy.py \
  --data "$DATA_DIR" \
  --max-seq-len "$MAX_SEQ_LEN" \
  --label-mode root_descendant \
  --also-freq-median \
  --out-dir "$OUT_DIR"

echo "Done: $(date)"
