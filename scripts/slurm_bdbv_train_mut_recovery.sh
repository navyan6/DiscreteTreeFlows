#!/bin/bash
#SBATCH --job-name=bdbv_mutrec
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/bdbv_mutrec_%j.log
#SBATCH --error=logs/bdbv_mutrec_%j.log
#
# BDBV L TreeSBM — covid_v5 mut-recipe fork.
#
# Examples:
#   sbatch --qos=mig-max scripts/slurm_bdbv_train_mut_recovery.sh
#   HOTSPOT=1 CKPT_DIR=checkpoints/bdbv_v1_lit_mutrec sbatch ...
#   DATA_ROOT=data/bdbv_temporal CKPT_DIR=checkpoints/bdbv_pan_v1_mutrec sbatch ...

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows

HOTSPOT="${HOTSPOT:-0}"
CKPT_DIR="${CKPT_DIR:-checkpoints/bdbv_v1_mutrec}"
DATA_ROOT="${DATA_ROOT:-data/bdbv_temporal}"
LAMBDA_BR="${LAMBDA_BR:-0}"
LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
LAMBDA_TOP="${LAMBDA_TOP:-0.1}"
MUT_HOTSPOT_MASK="${MUT_HOTSPOT_MASK:-results/bdbv_l_mask/mut_hotspot_mask_lit.pt}"
MUT_HOTSPOT_WEIGHT="${MUT_HOTSPOT_WEIGHT:-5.0}"
MUT_HOTSPOT_FORCE="${MUT_HOTSPOT_FORCE:-1}"

MAX_SEQ_LEN="${MAX_SEQ_LEN:-900}"
if [[ -f results/bdbv_l_conservation/window_config.json ]]; then
  MAX_SEQ_LEN=$($PYTHON -c "import json; print(json.load(open('results/bdbv_l_conservation/window_config.json'))['max_seq_len'])")
fi

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
  echo "WARN: empty val dir — using test band for early-stopping (sparse BDBV split)"
  VAL_DIR="$DATA_ROOT/test"
fi
echo "Start: $(date) ckpt=$CKPT_DIR data=$DATA_ROOT val=$VAL_DIR L=$MAX_SEQ_LEN hotspot=$HOTSPOT"

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
  --resume \
  "${EXTRA_ARGS[@]}"

echo "Done: $(date)"
