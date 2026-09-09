#!/bin/bash
#SBATCH --job-name=ab_oas_v2_shm
#SBATCH --partition=b200-mig90
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=8G
#SBATCH --time=72:00:00
#SBATCH --output=logs/ab_oas_v2_shm_%j.log
#SBATCH --error=logs/ab_oas_v2_shm_%j.log
#
# Track A v2: OAS TreeSBM with SHM Q0 site-rate prior (Recipe A).
# Does NOT overwrite ab_oas_1m_v1 or pathogen best.pt.
#
#   SKIP_PRECOMPUTE=1 REFRESH_CACHES=0 RESUME=1 \
#     sbatch --qos=mig-max scripts/slurm_ab_oas_v2_shm_train.sh

set -euo pipefail

export CKPT_DIR="${CKPT_DIR:-checkpoints/ab_oas_1m_v2_shm}"
export DATA="${DATA:-data/ab_clones_1m}"
export MUT_HOTSPOT_MASK="${MUT_HOTSPOT_MASK:-results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt}"
export MUT_HOTSPOT_WEIGHT="${MUT_HOTSPOT_WEIGHT:-8.0}"
export MUT_HOTSPOT_FORCE="${MUT_HOTSPOT_FORCE:-1}"
export SKIP_PRECOMPUTE="${SKIP_PRECOMPUTE:-1}"
export REFRESH_CACHES="${REFRESH_CACHES:-0}"
export RESUME="${RESUME:-1}"
export FITNESS_BETA="${FITNESS_BETA:-0}"
export SHM_SITE_BOOST="${SHM_SITE_BOOST:-2.0}"
export SHM_FWR_STAY="${SHM_FWR_STAY:-0.5}"
export SHM_USE_AID="${SHM_USE_AID:-1}"

# Refuse clobbering v1
if [[ "$CKPT_DIR" == "checkpoints/ab_oas_1m_v1" || "$CKPT_DIR" == "checkpoints/best" ]]; then
  echo "REFUSING overwrite of paper ckpt dir: $CKPT_DIR" >&2
  exit 2
fi

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows
mkdir -p logs "$CKPT_DIR"

# Reuse OAS train scaffolding via inline call (same data / entropy / CDR mask),
# but inject SHM prior flags into train.py.
MAX_SEQ_LEN="${MAX_SEQ_LEN:-160}"
LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
LAMBDA_BR="${LAMBDA_BR:-0}"
MUT_NORMALIZE="${MUT_NORMALIZE:-count}"
ENTROPY_ALPHA="${ENTROPY_ALPHA:-3.0}"
ENTROPY_ALPHA_CONS="${ENTROPY_ALPHA_CONS:-1.0}"
PSSM_GATE="${PSSM_GATE:-1}"
MUT_AA_EMB="${MUT_AA_EMB:-1}"
MUT_AA_EMB_DIM="${MUT_AA_EMB_DIM:-16}"

if [[ ! -d "$DATA/train" ]]; then
  echo "ERROR: missing $DATA/train" >&2
  exit 1
fi
if [[ ! -f "$MUT_HOTSPOT_MASK" ]]; then
  echo "ERROR: missing CDR mask $MUT_HOTSPOT_MASK" >&2
  exit 1
fi

EXTRA_ARGS=(
  --pssm-gate
  --mut-aa-emb --mut-aa-emb-dim "$MUT_AA_EMB_DIM"
  --entropy-weight-alpha-cons "$ENTROPY_ALPHA_CONS"
  --mut-hotspot-mask "$MUT_HOTSPOT_MASK"
  --mut-hotspot-weight "$MUT_HOTSPOT_WEIGHT"
  --mut-hotspot-force
  --shm-site-boost "$SHM_SITE_BOOST"
  --shm-fwr-stay "$SHM_FWR_STAY"
)
if [[ "$SHM_USE_AID" == "1" || "$SHM_USE_AID" == "true" ]]; then
  EXTRA_ARGS+=(--shm-use-aid)
fi
if [[ "$RESUME" == "1" || "$RESUME" == "true" ]]; then
  EXTRA_ARGS+=(--resume)
fi
# Explicit β=0 (do not pass --fitness-beta)
if [[ "$FITNESS_BETA" != "0" && "$FITNESS_BETA" != "0.0" ]]; then
  echo "REFUSING fitness β>0 for Ab SHM train (anti-SHM). Got FITNESS_BETA=$FITNESS_BETA" >&2
  exit 2
fi

echo "=== Ab OAS v2 SHM TreeSBM train ==="
echo "ckpt=$CKPT_DIR boost=$SHM_SITE_BOOST fwr_stay=$SHM_FWR_STAY use_aid=$SHM_USE_AID"
echo "hotspot_weight=$MUT_HOTSPOT_WEIGHT (v1 used 5; v2 default 8)"
ls -la "$MUT_HOTSPOT_MASK"

if [[ "$SKIP_PRECOMPUTE" != "1" ]]; then
  for split in train val; do
    $PYTHON -u scripts/precompute_plm.py --data "$DATA/$split"
    $PYTHON -u scripts/precompute_ref_rates.py --data "$DATA/$split" --max-seq-len "$MAX_SEQ_LEN"
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
