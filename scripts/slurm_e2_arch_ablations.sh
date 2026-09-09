#!/bin/bash
#SBATCH --job-name=e2_abl
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=72:00:00
#SBATCH --output=logs/e2_ablate_%j.log
#SBATCH --error=logs/e2_ablate_%j.log
#
# Appendix E.2 architecture ablations (COVID Brazil / covid_v5_mutrec recipe).
# Train-time: --ablate-mut-head / --ablate-stop-head (not gen-time flags on v5).
#
# ABLATE=no_mut_head|no_stop_head
# MODE=train (default) | enrich | coverage | baselines
#
#   ABLATE=no_mut_head MODE=train sbatch --qos=mig-max --time=72:00:00 \
#     scripts/slurm_e2_arch_ablations.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints benchmarks/results

ABLATE="${ABLATE:-no_mut_head}"
MODE="${MODE:-train}"
DATA="${DATA:-data/covid/test}"
TRAIN="${TRAIN:-data/covid/train}"
EVESCAPE=data/covid/evescape_spike_rbd.pt
MASK=results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt
MRS="${MRS:-0.5}"
NSTEPS="${NSTEPS:-100}"
MAX_TREES="${MAX_TREES:-20}"
K_MAX="${K_MAX:-100}"
K_STEP="${K_STEP:-10}"
MAX_ROOTS="${MAX_ROOTS:-5}"
BASE_MAX_ROOTS="${BASE_MAX_ROOTS:-20}"
BASE_K="${BASE_K:-20}"
BASE_M="${BASE_M:-20}"

case "$ABLATE" in
  no_mut_head)
    CKPT_DIR_DEFAULT=checkpoints/covid_e2_no_mut_head
    TRAIN_EXTRA=(--ablate-mut-head)
    ;;
  no_stop_head)
    CKPT_DIR_DEFAULT=checkpoints/covid_e2_no_stop_head
    TRAIN_EXTRA=(--ablate-stop-head)
    ;;
  *)
    echo "Unknown ABLATE=$ABLATE (no_mut_head|no_stop_head)"; exit 1
    ;;
esac

CKPT_DIR="${CKPT_DIR:-$CKPT_DIR_DEFAULT}"
CKPT="${CKPT:-${CKPT_DIR}/best.pt}"
TAG="e2_${ABLATE}"

resolve_ckpt() {
  local c="$1"
  if [[ -f "$c" ]]; then echo "$c"; return; fi
  local name; name=$(basename "$(dirname "$c")")
  local backup="$LABHOME/checkpoints_backup/${name}/best.pt"
  if [[ -f "$backup" ]]; then
    mkdir -p "$(dirname "$c")"
    cp -n "$backup" "$c" 2>/dev/null || cp "$backup" "$c"
    echo "$c"; return
  fi
  echo "$c"
}

if [[ "$MODE" == "train" ]]; then
  echo "E.2 train ablation: $ABLATE -> $CKPT_DIR  flags=${TRAIN_EXTRA[*]}"
  echo "Start: $(date)"
  export CKPT_DIR
  export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
  mkdir -p "$CKPT_DIR" logs
  LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
  LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
  LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
  LAMBDA_BR="${LAMBDA_BR:-0.1}"
  LAMBDA_TOP="${LAMBDA_TOP:-0.1}"
  MUT_NORMALIZE="${MUT_NORMALIZE:-count}"
  ENTROPY_ALPHA="${ENTROPY_ALPHA:-3.0}"
  ENTROPY_ALPHA_CONS="${ENTROPY_ALPHA_CONS:-1.0}"
  EPOCHS="${EPOCHS:-100}"
  PATIENCE="${PATIENCE:-50}"

  $PYTHON -u scripts/train.py \
      --data        data/covid/train \
      --val-data    data/covid/val \
      --test-data   data/covid/test \
      --max-seq-len 1280 \
      --epochs      "$EPOCHS" \
      --patience    "$PATIENCE" \
      --bridge-c    1.0 \
      --lambda-mut  "$LAMBDA_MUT" \
      --lambda-cons "$LAMBDA_CONS" \
      --lambda-br   "$LAMBDA_BR" \
      --lambda-top  "$LAMBDA_TOP" \
      --mut-normalize "$MUT_NORMALIZE" \
      --lambda-semi "$LAMBDA_SEMI" \
      --per-site-pos-emb \
      --use-site-entropy \
      --use-entropy-loss-weighting \
      --use-entropy-cons-weighting \
      --entropy-source empirical \
      --entropy-weight-alpha "$ENTROPY_ALPHA" \
      --entropy-weight-alpha-cons "$ENTROPY_ALPHA_CONS" \
      --entropy-weight-floor 1.0 \
      --pssm-gate \
      --mut-aa-emb \
      --ckpt-dir    "$CKPT_DIR" \
      --resume \
      "${TRAIN_EXTRA[@]}"
  echo "Done train: $(date) -> $CKPT_DIR/best.pt"
  exit 0
fi

CKPT_PATH=$(resolve_ckpt "$CKPT")
if [[ ! -f "$CKPT_PATH" ]]; then
  echo "Missing checkpoint $CKPT_PATH — train first (MODE=train)"; exit 1
fi

echo "Start: $(date) ablate=$ABLATE mode=$MODE ckpt=$CKPT_PATH"

if [[ "$MODE" == "enrich" ]]; then
  OUT="checkpoints/table_e2_${TAG}_mrs${MRS}.json"
  $PYTHON scripts/eval_evescape_enrichment.py \
    --checkpoint "$CKPT_PATH" \
    --data "$DATA" \
    --max-seq-len 1280 \
    --evescape "$EVESCAPE" \
    --lit-hotspot-mask "$MASK" \
    --mutation-rate-scale "$MRS" \
    --n-steps "$NSTEPS" \
    --max-trees "$MAX_TREES" \
    --out "$OUT"
  echo "Done: $(date) -> $OUT"
elif [[ "$MODE" == "coverage" ]]; then
  PARAMS=benchmarks/results/params_covid.json
  if [[ ! -f "$PARAMS" ]]; then
    $PYTHON benchmarks/fit_params.py --train-data "$TRAIN" --out "$PARAMS"
  fi
  OUT="benchmarks/results/coverage_curves_covid_N16_table_e2_${TAG}_K${K_MAX}.csv"
  CACHE_DIR="benchmarks/results/coverage_cache_covid_e2_${TAG}"
  $PYTHON benchmarks/coverage_curves.py \
    --test-data "$DATA" \
    --train-data "$TRAIN" \
    --params "$PARAMS" \
    --N 16 \
    --K-max "$K_MAX" --K-step "$K_STEP" \
    --max-roots "$MAX_ROOTS" \
    --methods treesbm \
    --eps-frac 0.02 \
    --e-list 0,1,2,3,5,8,10 \
    --checkpoint "$CKPT_PATH" \
    --n-steps 50 \
    --max-seq-len 1280 \
    --cache-dir "$CACHE_DIR" \
    --out "$OUT"
  echo "Done: $(date) -> $OUT"
elif [[ "$MODE" == "baselines" ]]; then
  PARAMS=benchmarks/results/params_covid.json
  if [[ ! -f "$PARAMS" ]]; then
    $PYTHON benchmarks/fit_params.py --train-data "$TRAIN" --out "$PARAMS"
  fi
  OUT="benchmarks/results/results_baselines_covid_e2_${TAG}.csv"
  $PYTHON benchmarks/run_table.py \
    --test-data "$DATA" \
    --train-data "$TRAIN" \
    --params "$PARAMS" \
    --N 16 --K "$BASE_K" --M "$BASE_M" --max-roots "$BASE_MAX_ROOTS" \
    --max-seq-len 1280 \
    --methods treesbm \
    --checkpoint "$CKPT_PATH" \
    --out "$OUT"
  echo "Done: $(date) -> $OUT"
else
  echo "Unknown MODE=$MODE (train|enrich|coverage|baselines)"; exit 1
fi
