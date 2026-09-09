#!/bin/bash
#SBATCH --job-name=hiv_train
#SBATCH --partition=b200-mig90
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/hiv_train_%j.log
#SBATCH --error=logs/hiv_train_%j.log
#
# HIV-1 Env TreeSBM train (lit-hotspot V1–V5 + entropy, BL head OFF).
#
# Usage:
#   DATA=data/hiv_temporal CKPT_DIR=checkpoints/hiv_temporal_v1 \
#     sbatch --dependency=afterok:<pipeline_job> scripts/slurm_hiv_train.sh
#   DATA=data/hiv_geo CKPT_DIR=checkpoints/hiv_geo_v1 \
#     sbatch --dependency=afterok:<pipeline_job> scripts/slurm_hiv_train.sh
#
# Locked settings:
#   LAMBDA_BR=0 (no BL prediction loss; TreeTime numdate bl.json still used for bridge times)
#   entropy ON (empirical)
#   mut hotspot = HXB2 V1–V5 mask
#   Env only, max_seq_len~900

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows
mkdir -p logs

DATA="${DATA:?set DATA=data/hiv_temporal or data/hiv_geo}"
CKPT_DIR="${CKPT_DIR:?set CKPT_DIR=checkpoints/hiv_temporal_v1 or checkpoints/hiv_geo_v1}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-900}"
LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
LAMBDA_BR="${LAMBDA_BR:-0}"
MUT_NORMALIZE="${MUT_NORMALIZE:-count}"
ENTROPY_ALPHA="${ENTROPY_ALPHA:-3.0}"
ENTROPY_ALPHA_CONS="${ENTROPY_ALPHA_CONS:-1.0}"
MUT_HOTSPOT_MASK="${MUT_HOTSPOT_MASK:-results/hiv_env_mask/mut_hotspot_mask_hiv_env_v1v5.pt}"
MUT_HOTSPOT_WEIGHT="${MUT_HOTSPOT_WEIGHT:-5.0}"
MUT_HOTSPOT_FORCE="${MUT_HOTSPOT_FORCE:-1}"
PSSM_GATE="${PSSM_GATE:-1}"
MUT_AA_EMB="${MUT_AA_EMB:-1}"
MUT_AA_EMB_DIM="${MUT_AA_EMB_DIM:-16}"
SKIP_PRECOMPUTE="${SKIP_PRECOMPUTE:-0}"
FITNESS_BETA="${FITNESS_BETA:-0}"
FITNESS_SCORE="${FITNESS_SCORE:-log_R0}"

mkdir -p "$CKPT_DIR"

if [[ ! -d "$DATA/train" ]]; then
  echo "ERROR: missing $DATA/train — run prepare_hiv_splits.py + pipeline first" >&2
  exit 1
fi

if [[ ! -f "$MUT_HOTSPOT_MASK" ]]; then
  echo "Building missing hotspot mask: $MUT_HOTSPOT_MASK"
  $PYTHON scripts/build_hiv_env_v1v5_hotspot_mask.py \
    --max-seq-len "$MAX_SEQ_LEN" \
    --out "$MUT_HOTSPOT_MASK"
fi

# Refuse non-HIV masks
case "$MUT_HOTSPOT_MASK" in
  *covid_mutfreq*|*pmc_lit*|*flu_mutfreq*|*nmicrobiol*|*ab_cdr*)
    echo "REFUSING non-HIV hotspot mask: $MUT_HOTSPOT_MASK" >&2
    exit 2
    ;;
esac

EXTRA_ARGS=()
if [[ "$PSSM_GATE" == "1" || "$PSSM_GATE" == "true" ]]; then
  EXTRA_ARGS+=(--pssm-gate)
fi
if [[ "$MUT_AA_EMB" == "1" || "$MUT_AA_EMB" == "true" ]]; then
  EXTRA_ARGS+=(--mut-aa-emb --mut-aa-emb-dim "$MUT_AA_EMB_DIM")
fi
if [[ "$FITNESS_BETA" != "0" && "$FITNESS_BETA" != "0.0" ]]; then
  EXTRA_ARGS+=(--fitness-beta "$FITNESS_BETA" --fitness-score "$FITNESS_SCORE")
fi
EXTRA_ARGS+=(--entropy-weight-alpha-cons "$ENTROPY_ALPHA_CONS")
EXTRA_ARGS+=(--mut-hotspot-mask "$MUT_HOTSPOT_MASK")
EXTRA_ARGS+=(--mut-hotspot-weight "$MUT_HOTSPOT_WEIGHT")
if [[ "$MUT_HOTSPOT_FORCE" == "1" || "$MUT_HOTSPOT_FORCE" == "true" ]]; then
  EXTRA_ARGS+=(--mut-hotspot-force)
fi

echo "=== HIV Env TreeSBM train ==="
echo "host=$(hostname) date=$(date -Is) job=${SLURM_JOB_ID:-local}"
echo "data=$DATA ckpt=$CKPT_DIR max_seq_len=$MAX_SEQ_LEN"
echo "lambda_br=$LAMBDA_BR (BL head loss disabled)"
echo "fitness_beta=$FITNESS_BETA fitness_score=$FITNESS_SCORE"
echo "entropy ON + hotspot=$MUT_HOTSPOT_MASK"
ls -la "$MUT_HOTSPOT_MASK"

for split in train val; do
  n=$($PYTHON -c "
from src.dataset import TreeDataset
ds = TreeDataset('${DATA}/$split', max_seq_len=int('${MAX_SEQ_LEN}'))
print(len(ds.groups))
" | tail -1)
  echo "TreeDataset $split: Found $n complete groups"
  if [[ "$n" -le 0 ]]; then
    echo "ERROR: 0 complete trees in $DATA/$split — pipeline not finished?" >&2
    exit 1
  fi
done

if [[ "$SKIP_PRECOMPUTE" != "1" ]]; then
  for split in train val test; do
    if [[ ! -d "$DATA/$split" ]]; then
      continue
    fi
    echo "=== precompute_plm $split ==="
    $PYTHON -u scripts/precompute_plm.py --data "$DATA/$split"
    echo "=== precompute_ref_rates $split ==="
    $PYTHON -u scripts/precompute_ref_rates.py \
      --data "$DATA/$split" \
      --max-seq-len "$MAX_SEQ_LEN"
  done
fi

echo "=== train.py ==="
$PYTHON -u scripts/train.py \
  --data "$DATA/train" \
  --val-data "$DATA/val" \
  --test-data "$DATA/test" \
  --max-seq-len "$MAX_SEQ_LEN" \
  --epochs 100 \
  --patience 50 \
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
  --ckpt-dir "$CKPT_DIR" \
  --resume \
  "${EXTRA_ARGS[@]}"

echo "Done: $(date -Is)"
ls -la "$CKPT_DIR/"
