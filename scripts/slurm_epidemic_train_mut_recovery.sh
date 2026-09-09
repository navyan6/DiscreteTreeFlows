#!/bin/bash
#SBATCH --job-name=epidemic_mutrec
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/epidemic_mutrec_%j.log
#SBATCH --error=logs/epidemic_mutrec_%j.log
#
# Generic mut-recipe train for epidemic/outbreak split dirs.
# Env: DATA_ROOT, CKPT_DIR, MAX_SEQ_LEN, optional HOTSPOT + MUT_HOTSPOT_MASK

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows

DATA_ROOT="${DATA_ROOT:?set DATA_ROOT}"
CKPT_DIR="${CKPT_DIR:?set CKPT_DIR}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"
HOTSPOT="${HOTSPOT:-0}"
LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
LAMBDA_BR="${LAMBDA_BR:-0}"
LAMBDA_TOP="${LAMBDA_TOP:-0.1}"
MUT_HOTSPOT_MASK="${MUT_HOTSPOT_MASK:-}"
MUT_HOTSPOT_WEIGHT="${MUT_HOTSPOT_WEIGHT:-5.0}"
MUT_HOTSPOT_FORCE="${MUT_HOTSPOT_FORCE:-1}"
INIT_CHECKPOINT="${INIT_CHECKPOINT:-}"
LR="${LR:-}"

EXTRA_ARGS=(--per-site-pos-emb --use-site-entropy --use-entropy-loss-weighting \
  --use-entropy-cons-weighting --entropy-source empirical --entropy-weight-alpha 3.0 \
  --entropy-weight-floor 1.0 --pssm-gate --mut-aa-emb --mut-aa-emb-dim 16)

if [[ "$HOTSPOT" == "1" || "$HOTSPOT" == "true" ]]; then
  if [[ -f "$MUT_HOTSPOT_MASK" ]]; then
    EXTRA_ARGS+=(--mut-hotspot-mask "$MUT_HOTSPOT_MASK")
    EXTRA_ARGS+=(--mut-hotspot-weight "$MUT_HOTSPOT_WEIGHT")
    if [[ "$MUT_HOTSPOT_FORCE" == "1" ]]; then
      EXTRA_ARGS+=(--mut-hotspot-force)
    fi
  else
    echo "WARNING: HOTSPOT=1 but missing $MUT_HOTSPOT_MASK" >&2
  fi
fi

mkdir -p logs "$CKPT_DIR"
VAL_DIR="$DATA_ROOT/val"
if [[ ! -d "$VAL_DIR" ]] || [[ -z "$(ls -A "$VAL_DIR"/*_group_*.fasta 2>/dev/null || true)" ]]; then
  echo "WARN: empty val — using test for early stopping"
  VAL_DIR="$DATA_ROOT/test"
fi

echo "Start: $(date) ckpt=$CKPT_DIR data=$DATA_ROOT val=$VAL_DIR L=$MAX_SEQ_LEN init=${INIT_CHECKPOINT:-none}"

TRAIN_ARGS=()
if [[ -n "$INIT_CHECKPOINT" ]]; then
  TRAIN_ARGS+=(--init-checkpoint "$INIT_CHECKPOINT")
elif [[ -f "$CKPT_DIR/best.pt" ]]; then
  TRAIN_ARGS+=(--resume)
fi

LR_ARGS=()
if [[ -n "$LR" ]]; then
  LR_ARGS+=(--lr "$LR")
fi

$PYTHON -u scripts/train.py \
  --data        "$DATA_ROOT/train" \
  --val-data    "$VAL_DIR" \
  --test-data   "$DATA_ROOT/test" \
  --max-seq-len "$MAX_SEQ_LEN" \
  --epochs      100 \
  --patience    50 \
  --bridge-c    1.0 \
  --lambda-mut  "$LAMBDA_MUT" \
  --lambda-cons "$LAMBDA_CONS" \
  --lambda-br   "$LAMBDA_BR" \
  --lambda-top  "$LAMBDA_TOP" \
  --mut-normalize count \
  --lambda-semi "$LAMBDA_SEMI" \
  --ckpt-dir    "$CKPT_DIR" \
  "${LR_ARGS[@]}" \
  "${TRAIN_ARGS[@]}" \
  "${EXTRA_ARGS[@]}"

echo "Done: $(date)"
