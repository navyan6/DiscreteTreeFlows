#!/bin/bash
#SBATCH --job-name=gen_leaf_pll
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/gen_leaf_pll_%j.log
#SBATCH --error=logs/gen_leaf_pll_%j.log
#
# Score PLL on all data/generated/** TreeSBM leaf sequences.
# Env: INPUT_DIR (default data/generated), BATCH_SIZE, LIMIT

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows
mkdir -p logs

INPUT_DIR="${INPUT_DIR:-data/generated}"
BATCH_SIZE="${BATCH_SIZE:-8}"
LIMIT_ARGS=()
if [[ -n "${LIMIT:-}" ]]; then
  LIMIT_ARGS+=(--limit "$LIMIT")
fi

echo "Start: $(date) input=$INPUT_DIR"
$PYTHON -u scripts/score_generated_leaf_pll.py \
  --input-dir "$INPUT_DIR" \
  --recursive \
  --r0-backend esm2 \
  --batch-size "$BATCH_SIZE" \
  "${LIMIT_ARGS[@]}"
echo "Done: $(date)"
