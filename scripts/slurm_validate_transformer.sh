#!/bin/bash
#SBATCH --job-name=transformer_val
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/transformer_val_%j.log
#SBATCH --error=logs/transformer_val_%j.log

set -e

CHECKPOINT="${1:?usage: sbatch scripts/slurm_validate_transformer.sh <checkpoint> <data_dir> [out_dir]}"
DATA_DIR="${2:?usage: sbatch scripts/slurm_validate_transformer.sh <checkpoint> <data_dir> [out_dir]}"
OUT_DIR="${3:-results/transformer_validation}"

mkdir -p logs "$OUT_DIR"

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows

echo "Start: $(date)"
echo "Checkpoint: $CHECKPOINT"
echo "Data: $DATA_DIR"

$PYTHON scripts/validate_transformer.py \
  --checkpoint "$CHECKPOINT" \
  --data "$DATA_DIR" \
  --out-dir "$OUT_DIR" \
  --device cuda \
  --max-seq-len 1280

echo "Done: $(date)"
