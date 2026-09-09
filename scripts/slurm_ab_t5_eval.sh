#!/bin/bash
#SBATCH --job-name=ab_t5_eval
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
# omit --mem: PARCC filter → 8GB/CPU on b200-mig45 (56GB/GPU)
#SBATCH --time=6:00:00
#SBATCH --output=logs/ab_t5_eval_%j.log
#SBATCH --error=logs/ab_t5_eval_%j.log
#
# Table 5 Ab baseline eval after preprocess (data/ab_clones/test).
# pLM needs GPU/ESM. Neutral SHM / AR report not_implemented unless flags set.
#
# Usage:
#   sbatch scripts/slurm_ab_t5_eval.sh
#   METHODS=plm_prior MAX_GROUPS=50 sbatch scripts/slurm_ab_t5_eval.sh
#   METHODS=neutral_shm,plm_prior,ar_tree_edit ALLOW_NEUTRAL_BD=1 sbatch ...

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs results/ab_t5

DATA="${DATA:-data/ab_clones/test}"
METHODS="${METHODS:-plm_prior}"
OUT="${OUT:-results/ab_t5/eval_${METHODS//,/_}.json}"
MAX_GROUPS="${MAX_GROUPS:-50}"
K="${K:-10}"
ALLOW_NEUTRAL_BD="${ALLOW_NEUTRAL_BD:-0}"

echo "Start: $(date)"
echo "DATA=$DATA METHODS=$METHODS OUT=$OUT"

if [[ ! -d "$DATA" ]]; then
  echo "ERROR: missing $DATA — run slurm_ab_t5_preprocess.sh first"
  exit 1
fi

EXTRA=()
if [[ "$ALLOW_NEUTRAL_BD" == "1" ]]; then
  EXTRA+=(--allow-neutral-bd-fallback)
fi

$PYTHON scripts/eval_ab_maturation.py \
  --data "$DATA" \
  --methods "$METHODS" \
  --out "$OUT" \
  --max-groups "$MAX_GROUPS" \
  --K "$K" \
  "${EXTRA[@]}"

echo "Done: $(date) → $OUT"
