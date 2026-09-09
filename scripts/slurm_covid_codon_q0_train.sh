#!/bin/bash
#SBATCH --job-name=covid_codon_q0
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=72:00:00
#SBATCH --output=logs/covid_codon_q0_%j.log
#SBATCH --error=logs/covid_codon_q0_%j.log
#
# COVID codon-GY94 Q0 + learned residual c_θ, β=0.
# New dir only. Refuses paper / mutlin / OAS ckpt paths.

set -euo pipefail

export CKPT_DIR="${CKPT_DIR:-checkpoints/covid_v8_codon_q0}"
export DATA="${DATA:-data/covid}"
export R0_BACKEND="${R0_BACKEND:-codon_gy94}"
export FITNESS_BETA="${FITNESS_BETA:-0}"
export SKIP_PRECOMPUTE="${SKIP_PRECOMPUTE:-0}"
export RESUME="${RESUME:-1}"

for keep in \
  checkpoints/covid_v5_mutrec checkpoints/h1n1_v2_lit_hotspot \
  checkpoints/h3n2_v3_lit_hotspot checkpoints/hiv_geo_v1 \
  checkpoints/ab_oas_1m_v1 checkpoints/ab_oas_1m_v2_shm \
  checkpoints/ab_oas_1m_v3_thrifty checkpoints/best \
  checkpoints/covid_v6_mutlin_b0 checkpoints/covid_v6_mutlin_b025
 do
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
mkdir -p logs "$CKPT_DIR"

MAX_SEQ_LEN="${MAX_SEQ_LEN:-1280}"
LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
LAMBDA_BR="${LAMBDA_BR:-0.1}"
MUT_NORMALIZE="${MUT_NORMALIZE:-count}"

if [[ "$SKIP_PRECOMPUTE" != "1" ]]; then
  for split in train val test; do
    $PYTHON -u scripts/precompute_ref_rates.py \
      --data "$DATA/$split" --max-seq-len "$MAX_SEQ_LEN" --r0-backend "$R0_BACKEND"
  done
fi

EXTRA=(--r0-backend "$R0_BACKEND" --per-site-pos-emb --use-site-entropy
       --use-entropy-loss-weighting --use-entropy-cons-weighting
       --entropy-source empirical --entropy-weight-alpha 3.0 --entropy-weight-floor 1.0)
if [[ "$RESUME" == "1" ]]; then EXTRA+=(--resume); fi
if [[ "$FITNESS_BETA" != "0" && "$FITNESS_BETA" != "0.0" ]]; then
  echo "REFUSING FITNESS_BETA=$FITNESS_BETA (codon Q0 recipe is β=0)" >&2
  exit 2
fi

$PYTHON -u scripts/train.py \
  --data "$DATA/train" --val-data "$DATA/val" --test-data "$DATA/test" \
  --max-seq-len "$MAX_SEQ_LEN" --epochs 100 --patience 50 --lr 1e-4 \
  --bridge-c 1.0 --lambda-mut "$LAMBDA_MUT" --lambda-cons "$LAMBDA_CONS" \
  --lambda-top 0.1 --lambda-br "$LAMBDA_BR" --lambda-semi "$LAMBDA_SEMI" \
  --mut-normalize "$MUT_NORMALIZE" --ckpt-dir "$CKPT_DIR" \
  "${EXTRA[@]}"

echo "Done $(date -Is)"
ls -la "$CKPT_DIR/"
