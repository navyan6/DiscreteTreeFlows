#!/bin/bash
#SBATCH --job-name=eval_bdbv
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/eval_bdbv_%j.log
#SBATCH --error=logs/eval_bdbv_%j.log
#
# CKPT=checkpoints/bdbv_v1_mutrec sbatch --qos=mig-max scripts/slurm_eval_bdbv.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

DATA_ROOT="${DATA_ROOT:-data/bdbv_temporal}"
CKPT="${CKPT:-checkpoints/bdbv_v1_mutrec/best.pt}"
CKPT_PATH="${CKPT}/best.pt"
[[ -f "$CKPT" ]] && CKPT_PATH="$CKPT"
CKPT_NAME=$(basename "$(dirname "$CKPT_PATH")")
MAX_SEQ_LEN="${MAX_SEQ_LEN:-900}"
if [[ -f results/bdbv_l_conservation/window_config.json ]]; then
  MAX_SEQ_LEN=$($PYTHON -c "import json; print(json.load(open('results/bdbv_l_conservation/window_config.json'))['max_seq_len'])")
fi

EVAL_MRS="${EVAL_MRS:-0.3 0.5 1.0}"
EVAL_MAX_TREES="${EVAL_MAX_TREES:-20}"
NSTEPS="${NSTEPS:-100}"
LIT_MASK="${LIT_MASK:-results/bdbv_l_mask/mut_hotspot_mask_lit.pt}"
LIT_ARGS=()
if [[ -f "$LIT_MASK" ]]; then
  LIT_ARGS=(--lit-hotspot-mask "$LIT_MASK")
fi

for MRS in $EVAL_MRS; do
  OUT="checkpoints/eval_enrichment_${CKPT_NAME}_mrs${MRS}.json"
  echo "enrichment $CKPT_NAME mrs=$MRS -> $OUT"
  $PYTHON scripts/eval_evescape_enrichment.py \
    --checkpoint "$CKPT_PATH" \
    --data "$DATA_ROOT/test" \
    --max-seq-len "$MAX_SEQ_LEN" \
    --mutation-rate-scale "$MRS" \
    --n-steps "$NSTEPS" \
    --max-trees "$EVAL_MAX_TREES" \
    --out "$OUT" \
    "${LIT_ARGS[@]}"
done

VAR_OUT="checkpoints/eval_bdbv_2026_variants_${CKPT_NAME}.json"
$PYTHON scripts/eval_bdbv_outbreak_variants.py \
  --checkpoint "$CKPT_PATH" \
  --data "$DATA_ROOT/test" \
  --max-seq-len "$MAX_SEQ_LEN" \
  --mutation-rate-scale 0.5 \
  --n-steps "$NSTEPS" \
  --max-trees "$EVAL_MAX_TREES" \
  --out "$VAR_OUT"

echo "Done: $(date)"
