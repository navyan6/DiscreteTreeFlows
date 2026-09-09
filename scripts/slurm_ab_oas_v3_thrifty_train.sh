#!/bin/bash
#SBATCH --job-name=ab_oas_v3_thr
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=72:00:00
#SBATCH --output=logs/ab_oas_v3_thrifty_train_%j.log
#SBATCH --error=logs/ab_oas_v3_thrifty_train_%j.log
#
# Recipe B: OAS TreeSBM with Thrifty AA Q0, β=0, no Recipe A SHM stay prior.
# Does NOT overwrite ab_oas_1m_v1, ab_oas_1m_v2_shm, or pathogen best.pt.
#
#   SKIP_PRECOMPUTE=1 sbatch --dependency=afterok:<pre_job> \
#     scripts/slurm_ab_oas_v3_thrifty_train.sh

set -euo pipefail

export CKPT_DIR="${CKPT_DIR:-checkpoints/ab_oas_1m_v3_thrifty}"
export DATA="${DATA:-data/ab_clones_1m}"
export MUT_HOTSPOT_MASK="${MUT_HOTSPOT_MASK:-results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt}"
export MUT_HOTSPOT_WEIGHT="${MUT_HOTSPOT_WEIGHT:-5.0}"
export MUT_HOTSPOT_FORCE="${MUT_HOTSPOT_FORCE:-1}"
export SKIP_PRECOMPUTE="${SKIP_PRECOMPUTE:-1}"
export RESUME="${RESUME:-1}"
export FITNESS_BETA="${FITNESS_BETA:-0}"
# Antibody Recipe B only. Viral / OAS v1/v2 jobs must keep default esm2.
export R0_BACKEND="${R0_BACKEND:-thrifty_aa}"

for keep in checkpoints/ab_oas_1m_v1 checkpoints/ab_oas_1m_v2_shm checkpoints/best; do
  if [[ "$CKPT_DIR" == "$keep" ]]; then
    echo "REFUSING overwrite of locked ckpt dir: $CKPT_DIR" >&2
    exit 2
  fi
done

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows
export PYTHONPATH="$PWD:${PWD}/antibody_benchmark/data/raw/repos/netam:${PYTHONPATH:-}"
mkdir -p logs "$CKPT_DIR"

MAX_SEQ_LEN="${MAX_SEQ_LEN:-160}"
LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
LAMBDA_BR="${LAMBDA_BR:-0}"
MUT_NORMALIZE="${MUT_NORMALIZE:-count}"
ENTROPY_ALPHA="${ENTROPY_ALPHA:-3.0}"
ENTROPY_ALPHA_CONS="${ENTROPY_ALPHA_CONS:-1.0}"

if [[ ! -d "$DATA/train" ]]; then
  echo "ERROR: missing $DATA/train" >&2
  exit 1
fi
if [[ ! -f "$MUT_HOTSPOT_MASK" ]]; then
  echo "ERROR: missing CDR mask $MUT_HOTSPOT_MASK" >&2
  exit 1
fi

n_thr=$(ls -1 "$DATA/train"/group_*_ref_rates_thrifty.pt 2>/dev/null | wc -l | tr -d ' ')
if [[ "$n_thr" -le 0 ]]; then
  echo "ERROR: no group_*_ref_rates_thrifty.pt under $DATA/train — run precompute first" >&2
  exit 1
fi

EXTRA_ARGS=(
  --pssm-gate
  --mut-aa-emb --mut-aa-emb-dim 16
  --entropy-weight-alpha-cons "$ENTROPY_ALPHA_CONS"
  --mut-hotspot-mask "$MUT_HOTSPOT_MASK"
  --mut-hotspot-weight "$MUT_HOTSPOT_WEIGHT"
  --mut-hotspot-force
  --r0-backend "$R0_BACKEND"
)
if [[ "$RESUME" == "1" || "$RESUME" == "true" ]]; then
  EXTRA_ARGS+=(--resume)
fi
if [[ "$FITNESS_BETA" != "0" && "$FITNESS_BETA" != "0.0" ]]; then
  echo "REFUSING fitness β>0 for Recipe B (anti-SHM). Got FITNESS_BETA=$FITNESS_BETA" >&2
  exit 2
fi

echo "=== Ab OAS v3 Thrifty Q0 TreeSBM train ==="
echo "ckpt=$CKPT_DIR r0=$R0_BACKEND thrifty_caches=$n_thr"
ls -la "$MUT_HOTSPOT_MASK"

if [[ "$SKIP_PRECOMPUTE" != "1" ]]; then
  for split in train val; do
    $PYTHON -u scripts/precompute_ref_rates.py \
      --data "$DATA/$split" --max-seq-len "$MAX_SEQ_LEN" --r0-backend "$R0_BACKEND"
  done
fi

$PYTHON -u scripts/train.py \
  --data "$DATA/train" \
  --val-data "$DATA/val" \
  --test-data "$DATA/test" \
  --max-seq-len "$MAX_SEQ_LEN" \
  --epochs 150 \
  --patience 40 \
  --lr 1e-4 \
  --bridge-c 1.0 \
  --lambda-mut "$LAMBDA_MUT" \
  --lambda-cons "$LAMBDA_CONS" \
  --lambda-top 0.5 \
  --lambda-br "$LAMBDA_BR" \
  --lambda-semi "$LAMBDA_SEMI" \
  --mut-normalize "$MUT_NORMALIZE" \
  --per-site-pos-emb \
  --use-site-entropy \
  --use-entropy-loss-weighting \
  --use-entropy-cons-weighting \
  --entropy-source empirical \
  --entropy-weight-alpha "$ENTROPY_ALPHA" \
  --entropy-weight-floor 1.0 \
  --n-t-samples 4 \
  --ckpt-dir "$CKPT_DIR" \
  "${EXTRA_ARGS[@]}"

echo "Done: $(date -Is)"
ls -la "$CKPT_DIR/"
