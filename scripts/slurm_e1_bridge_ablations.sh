#!/bin/bash
#SBATCH --job-name=e1_abl
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/e1_ablate_%j.log
#SBATCH --error=logs/e1_ablate_%j.log
#
# Appendix E.1 bridge-matching train ablations (COVID Brazil / covid_v5_mutrec recipe).
# These are NOT gen-time --ablate-* on covid_v5_mutrec — they change the train target.
#
# ABLATE=terminal_only|no_doob|no_terminal
# MODE=train (default) | enrich | coverage | baselines
#
#   ABLATE=terminal_only MODE=train sbatch --qos=mig-max --time=24:00:00 \
#     scripts/slurm_e1_bridge_ablations.sh
#   ABLATE=no_doob MODE=enrich CKPT=checkpoints/covid_e1_no_doob/best.pt \
#     sbatch --qos=mig-max scripts/slurm_e1_bridge_ablations.sh
#
# After train jobs finish, re-submit with MODE=enrich|coverage and matching CKPT.

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints benchmarks/results

ABLATE="${ABLATE:-terminal_only}"
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
  terminal_only)
    CKPT_DIR_DEFAULT=checkpoints/covid_e1_terminal_only
    TRAIN_EXTRA=(--ablate-terminal-only)
    ;;
  no_doob)
    CKPT_DIR_DEFAULT=checkpoints/covid_e1_no_doob
    TRAIN_EXTRA=(--ablate-doob)
    ;;
  no_terminal)
    CKPT_DIR_DEFAULT=checkpoints/covid_e1_no_terminal
    TRAIN_EXTRA=(--ablate-terminal-consistency)
    ;;
  *)
    echo "Unknown ABLATE=$ABLATE (terminal_only|no_doob|no_terminal)"; exit 1
    ;;
esac

CKPT_DIR="${CKPT_DIR:-$CKPT_DIR_DEFAULT}"
CKPT="${CKPT:-${CKPT_DIR}/best.pt}"
TAG="e1_${ABLATE}"

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
  echo "E.1 train ablation: $ABLATE -> $CKPT_DIR  flags=${TRAIN_EXTRA[*]}"
  echo "Start: $(date)"
  # Same covid_v5_mutrec recipe; do NOT resume from full v5 (different loss).
  # --resume only applies if this ablation dir already has best.pt (time-limit restart).
  export CKPT_DIR
  # Inject ablation flags into EXTRA_ARGS used by the train recipe.
  # slurm_covid_train_mut_recovery.sh appends EXTRA_ARGS; we re-exec via env.
  # Prefer direct train.py call matching the v5 recipe so flags are explicit.
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
  OUT="checkpoints/table_e1_${TAG}_mrs${MRS}.json"
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
  OUT="benchmarks/results/coverage_curves_covid_N16_table_e1_${TAG}_K${K_MAX}.csv"
  CACHE_DIR="benchmarks/results/coverage_cache_covid_e1_${TAG}"
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
  OUT="benchmarks/results/results_baselines_covid_e1_${TAG}.csv"
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
