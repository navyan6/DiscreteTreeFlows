#!/bin/bash
#SBATCH --job-name=bdbv_cov
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/bdbv_coverage_%j.log
#SBATCH --error=logs/bdbv_coverage_%j.log
#
# CHECKPOINT=checkpoints/bdbv_v1_mutrec/best.pt sbatch --qos=mig-max scripts/slurm_bdbv_coverage.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results

DATA_ROOT="${DATA_ROOT:-data/bdbv_temporal}"
CHECKPOINT="${CHECKPOINT:-checkpoints/bdbv_v1_mutrec/best.pt}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-900}"
if [[ -f results/bdbv_l_conservation/window_config.json ]]; then
  MAX_SEQ_LEN=$($PYTHON -c "import json; print(json.load(open('results/bdbv_l_conservation/window_config.json'))['max_seq_len'])")
fi

TEST_DATA="$DATA_ROOT/test"
TRAIN_DATA="$DATA_ROOT/train"
PARAMS="${PARAMS:-benchmarks/results/params_bdbv.json}"
OUT="${OUT:-benchmarks/results/coverage_curves_bdbv_N16_eabs.csv}"
CACHE_DIR="${CACHE_DIR:-benchmarks/results/coverage_cache_bdbv_eabs}"
METHODS="${METHODS:-neutral_bd plm_prior artreeformer_adapted treesbm}"

if [[ ! -f "$PARAMS" ]]; then
  $PYTHON benchmarks/fit_params.py --train-data "$TRAIN_DATA" --out "$PARAMS"
fi

$PYTHON benchmarks/coverage_curves.py \
  --test-data "$TEST_DATA" \
  --train-data "$TRAIN_DATA" \
  --params "$PARAMS" \
  --N 16 \
  --K-max 100 --K-step 10 \
  --max-roots 5 \
  --methods $METHODS \
  --eps-frac 0.02 \
  --e-list 0,1,2,3,5,8,10 \
  --checkpoint "$CHECKPOINT" \
  --n-steps 50 \
  --max-seq-len "$MAX_SEQ_LEN" \
  --out "$OUT" \
  --cache-dir "$CACHE_DIR"

echo "Done: $(date) -> $OUT"
