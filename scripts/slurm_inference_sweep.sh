#!/bin/bash
#SBATCH --job-name=inf_sweep
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/inf_sweep_%j.log
#SBATCH --error=logs/inf_sweep_%j.log
#
# High-mrs inference sweep (mrs 0.3/0.5/1.0) + optional site-propensity sampling.
# Skips finished tags. Coordinates with Agent G: default ckpt=v3; override for v4/v5.
#
#   sbatch --qos=mig-max scripts/slurm_inf_sweep.sh
#   CKPT=checkpoints/covid_v4_mutrec/best.pt sbatch --qos=mig-max scripts/slurm_inf_sweep.sh
#   SITE_SOFTMAX=1 CKPT=checkpoints/covid_v3_cons/best.pt \
#     sbatch --qos=mig-max scripts/slurm_inf_sweep.sh
#   EVAL_MRS="0.3 0.5 1.0" NSTEPS_LIST="100" MAX_TREES=20 \
#     sbatch --qos=mig-max scripts/slurm_inf_sweep.sh

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

CKPT="${CKPT:-checkpoints/covid_v3_cons/best.pt}"
if [[ ! -f "$CKPT" && -f "$LABHOME/checkpoints_backup/covid_v3_cons/best.pt" ]]; then
    mkdir -p checkpoints/covid_v3_cons
    cp -n "$LABHOME/checkpoints_backup/covid_v3_cons/best.pt" checkpoints/covid_v3_cons/best.pt
    CKPT=checkpoints/covid_v3_cons/best.pt
fi
if [[ ! -f "$CKPT" ]]; then
    echo "ERROR: checkpoint not found: $CKPT"
    exit 1
fi

MRS_LIST="${MRS_LIST:-${EVAL_MRS:-0.3 0.5 1.0}}"
NSTEPS_LIST="${NSTEPS_LIST:-50 100 150}"
MAX_TREES="${MAX_TREES:-10}"
SKIP_EXISTING="${SKIP_EXISTING:-1}"
SITE_SOFTMAX="${SITE_SOFTMAX:-0}"
SITE_TEMPERATURE="${SITE_TEMPERATURE:-1.0}"

SUFFIX=""
SOFTMAX_ARGS=()
if [[ "$SITE_SOFTMAX" == "1" || "$SITE_SOFTMAX" == "true" ]]; then
    SUFFIX="_sitesm"
    SOFTMAX_ARGS+=(--site-softmax-sample --site-temperature "$SITE_TEMPERATURE")
fi

OUT_DIR="checkpoints/inference_sweep_$(basename "$(dirname "$CKPT")")${SUFFIX}"
mkdir -p "$OUT_DIR"
echo "Start: $(date)  ckpt=$CKPT  out=$OUT_DIR"
echo "MRS_LIST=$MRS_LIST  NSTEPS_LIST=$NSTEPS_LIST  MAX_TREES=$MAX_TREES"
echo "SKIP_EXISTING=$SKIP_EXISTING  SITE_SOFTMAX=$SITE_SOFTMAX"

for MRS in $MRS_LIST; do
  for NSTEPS in $NSTEPS_LIST; do
    TAG="mrs${MRS}_steps${NSTEPS}"
    OUT="$OUT_DIR/eval_${TAG}.json"
    if [[ "$SKIP_EXISTING" == "1" && -f "$OUT" ]]; then
      echo "=== SKIP $TAG (exists: $OUT) ==="
      continue
    fi
    echo "=== $TAG ==="
    $PYTHON scripts/eval_evescape_enrichment.py \
        --checkpoint "$CKPT" \
        --data data/covid/test \
        --max-seq-len 1280 \
        --evescape data/covid/evescape_spike_rbd.pt \
        --mutation-rate-scale "$MRS" \
        --n-steps "$NSTEPS" \
        --max-trees "$MAX_TREES" \
        --out "$OUT" \
        "${SOFTMAX_ARGS[@]}"
  done
done

echo "Done: $(date)"
echo "Pick config with best mut_recovery subject to cons_retention >= 0.98"
echo "Results under $OUT_DIR"
