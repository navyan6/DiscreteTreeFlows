#!/bin/bash
#SBATCH --job-name=transformer_val
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/transformer_val_%j.log
#SBATCH --error=logs/transformer_val_%j.log
#
# Frozen GraphTF linear-probe validation (Appendix Table J.1).
# Default QOS: mig-max (avoids lab MaxGRESPerAccount on qos=mig).
# Override: sbatch --qos=<qos> scripts/slurm_validate_transformer.sh ...
#   or:     QOS=mig-max sbatch --qos="${QOS}" scripts/slurm_validate_transformer.sh ...
#
# Usage:
#   sbatch scripts/slurm_validate_transformer.sh <checkpoint> <data_dir> [out_dir] [max_seq_len] [max_trees]
#   sbatch scripts/slurm_validate_transformer.sh checkpoints/h3n2_v2/best.pt data/h3n2/train results/transformer_val_h3n2_j1 566
#   sbatch scripts/slurm_validate_transformer.sh checkpoints/covid_v3_cons/best.pt data/covid/train results/transformer_val_covid 1280

set -e

CHECKPOINT="${1:?usage: sbatch scripts/slurm_validate_transformer.sh <checkpoint> <data_dir> [out_dir] [max_seq_len] [max_trees]}"
DATA_DIR="${2:?usage: sbatch scripts/slurm_validate_transformer.sh <checkpoint> <data_dir> [out_dir] [max_seq_len] [max_trees]}"
OUT_DIR="${3:-results/transformer_validation}"
MAX_SEQ_LEN="${4:-566}"
MAX_TREES="${5:-}"

mkdir -p logs "$OUT_DIR"

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_home}"

cd ~/DiscreteTreeFlows

echo "Start: $(date)"
echo "Checkpoint: $CHECKPOINT"
echo "Data: $DATA_DIR"
echo "Out: $OUT_DIR  max_seq_len=$MAX_SEQ_LEN  max_trees=${MAX_TREES:-all}"

EXTRA=()
if [[ -n "$MAX_TREES" ]]; then
  EXTRA+=(--max-trees "$MAX_TREES")
fi

$PYTHON scripts/validate_transformer.py \
  --checkpoint "$CHECKPOINT" \
  --data "$DATA_DIR" \
  --out-dir "$OUT_DIR" \
  --device cuda \
  --max-seq-len "$MAX_SEQ_LEN" \
  "${EXTRA[@]}"

echo "Done: $(date)"
echo "Primary J.1 table: $OUT_DIR/probe_table_j1_primary.csv"
