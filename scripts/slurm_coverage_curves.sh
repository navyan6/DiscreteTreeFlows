#!/bin/bash
#SBATCH --job-name=cov_curves
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/coverage_curves_%j.log
#SBATCH --error=logs/coverage_curves_%j.log
#
# Table 3 coverage curves (NeutralBD / PLMPrior / artreeformer_adapted / TreeSBM).
# Locked: site_recall, eps=0.02, 5 H3N2 roots, K=10..100 step 10.
# Prefer a small smoke first before the full 4-method sweep.
#
# Smoke (NeutralBD only):
#   METHODS="neutral_bd" K_MAX=20 MAX_ROOTS=2 N_LIST="16" NO_ESM=1 \
#     OUT=benchmarks/results/coverage_curves_smoke.csv \
#     sbatch --qos=mig-max --time=1:00:00 scripts/slurm_coverage_curves.sh
#
# Full locked submit (5 roots × 4 methods × K_max=100) WITH tree cache:
#   METHODS="neutral_bd plm_prior artreeformer_adapted treesbm" \
#     K_MAX=100 K_STEP=10 MAX_ROOTS=5 N_LIST="16" \
#     CHECKPOINT=checkpoints/h3n2_v3_lit_hotspot/best.pt \
#     E_LIST="0,1,2,3,5,8,10" \
#     CACHE_DIR=benchmarks/results/coverage_cache_h3n2_N16 \
#     OUT=benchmarks/results/coverage_curves_h3n2_N16_eabs.csv \
#     sbatch --qos=mig-max scripts/slurm_coverage_curves.sh
#
# Rescore-only (after CACHE_DIR exists; no TreeSBM regen):
#   RESCORE_FROM=benchmarks/results/coverage_cache_h3n2_N16 \
#     METHODS="neutral_bd plm_prior artreeformer_adapted treesbm" \
#     E_LIST="0,1,2,3,5" \
#     OUT=benchmarks/results/coverage_curves_h3n2_N16_eabs_rescore.csv \
#     sbatch --qos=mig-max --time=1:00:00 scripts/slurm_coverage_curves.sh
#
# Prereqs on Betty:
#   - data/h3n2/{train,test}/group_*_rooted.nwk + *_anc_aa.fasta
#   - benchmarks/results/params.json  (or fit_params.py runs below)
#   - for artreeformer_adapted: benchmarks/external_pools/sampled/artreeformer_N16.nwk
#   - for treesbm: CHECKPOINT (default h3n2_v3_lit_hotspot/best.pt)
# NOTE: job 7439926 did NOT save trees — cannot rescore that run; need CACHE_DIR.

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results/plots

TEST_DATA="${TEST_DATA:-data/h3n2/test}"
TRAIN_DATA="${TRAIN_DATA:-data/h3n2/train}"
PARAMS="${PARAMS:-benchmarks/results/params.json}"
N_LIST="${N_LIST:-16}"
K_MAX="${K_MAX:-100}"
K_STEP="${K_STEP:-10}"
MAX_ROOTS="${MAX_ROOTS:-5}"
METHODS="${METHODS:-neutral_bd plm_prior artreeformer_adapted treesbm}"
EPS_FRAC="${EPS_FRAC:-0.02}"
E_LIST="${E_LIST:-0,1,2,3,5,8,10}"
CACHE_DIR="${CACHE_DIR:-}"
RESCORE_FROM="${RESCORE_FROM:-}"
CHECKPOINT="${CHECKPOINT:-checkpoints/h3n2_v3_lit_hotspot/best.pt}"
N_STEPS="${N_STEPS:-50}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"
OUT="${OUT:-benchmarks/results/coverage_curves_h3n2_N16_eabs.csv}"
NO_ESM="${NO_ESM:-0}"
MIN_CLADE_SIZE="${MIN_CLADE_SIZE:-2}"
MIN_SHARED_MUTS="${MIN_SHARED_MUTS:-1}"

echo "Start: $(date)"
echo "test=$TEST_DATA train=$TRAIN_DATA N=($N_LIST) K_max=$K_MAX step=$K_STEP"
echo "max_roots=$MAX_ROOTS methods=($METHODS) eps=$EPS_FRAC e_list=$E_LIST"
echo "checkpoint=$CHECKPOINT out=$OUT"
echo "cache_dir=${CACHE_DIR:-none} rescore_from=${RESCORE_FROM:-none}"

if [ ! -f "$PARAMS" ]; then
  echo "Fitting BD params on $TRAIN_DATA -> $PARAMS"
  $PYTHON benchmarks/fit_params.py --train-data "$TRAIN_DATA" --out "$PARAMS"
fi

EXTRA=()
if [ "$NO_ESM" = "1" ]; then
  EXTRA+=(--no-esm)
fi
if [ -n "$CACHE_DIR" ]; then
  EXTRA+=(--cache-dir "$CACHE_DIR")
fi
if [ -n "$RESCORE_FROM" ]; then
  EXTRA+=(--rescore-from "$RESCORE_FROM")
fi

# shellcheck disable=SC2086
$PYTHON benchmarks/coverage_curves.py \
  --test-data "$TEST_DATA" \
  --train-data "$TRAIN_DATA" \
  --params "$PARAMS" \
  --N $N_LIST \
  --K-max "$K_MAX" \
  --K-step "$K_STEP" \
  --max-roots "$MAX_ROOTS" \
  --methods $METHODS \
  --eps-frac "$EPS_FRAC" \
  --e-list "$E_LIST" \
  --min-clade-size "$MIN_CLADE_SIZE" \
  --min-shared-muts "$MIN_SHARED_MUTS" \
  --checkpoint "$CHECKPOINT" \
  --n-steps "$N_STEPS" \
  --max-seq-len "$MAX_SEQ_LEN" \
  --out "$OUT" \
  "${EXTRA[@]}"

if command -v Rscript >/dev/null 2>&1; then
  Rscript benchmarks/plot_coverage_curves.R "$OUT" \
    "benchmarks/results/plots/$(basename "${OUT%.csv}").pdf" || \
    echo "R plot skipped (Rscript failed)"
else
  echo "Rscript not found — CSV ready for offline plot"
fi

echo "Done: $(date)"
echo "CSV: $OUT"
