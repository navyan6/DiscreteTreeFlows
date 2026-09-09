#!/bin/bash
#SBATCH --job-name=c3_horizon
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=18:00:00
#SBATCH --output=logs/table_c3_horizon_%j.log
#SBATCH --error=logs/table_c3_horizon_%j.log
#
# Table C.3 — genetic H tercile stratification (true conditioning).
#
#   VIRUS=covid sbatch --qos=mig-max scripts/slurm_table_c3_horizon.sh
#   VIRUS=h1n1  sbatch --qos=mig-max scripts/slurm_table_c3_horizon.sh
#   VIRUS=h3n2  sbatch --qos=mig-max scripts/slurm_table_c3_horizon.sh
#   VIRUS=hiv_geo sbatch --qos=mig-max scripts/slurm_table_c3_horizon.sh
#   VIRUS=hiv_temporal sbatch --qos=mig-max scripts/slurm_table_c3_horizon.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results/tables

VIRUS="${VIRUS:-covid}"
K_LIST="${K_LIST:-10 100}"
MAX_ROOTS="${MAX_ROOTS:-5}"
N="${N:-16}"
E_LIST="${E_LIST:-1,2,3,5}"

# Defaults match table_c3_horizon.VIRUS_DEFAULTS / Table 5+4 roots.
case "$VIRUS" in
  covid)
    PARAMS=benchmarks/results/params_covid.json
    TRAIN_DATA=data/covid/train
    CHECKPOINT="${CHECKPOINT:-checkpoints/covid_v5_mutrec/best.pt}"
    MAX_SEQ_LEN=1280
    PIN_GROUPS="${PIN_GROUPS:-4,5,7,9,10}"
    METHODS="${METHODS:-neutral_bd artreeformer_adapted treesbm}"
    ;;
  h1n1)
    PARAMS=benchmarks/results/params_h1n1.json
    TRAIN_DATA=data/h1n1/train
    CHECKPOINT="${CHECKPOINT:-checkpoints/h1n1_v2_lit_hotspot/best.pt}"
    MAX_SEQ_LEN=566
    PIN_GROUPS="${PIN_GROUPS:-1,4,5,6,7}"
    METHODS="${METHODS:-neutral_bd artreeformer_adapted treesbm}"
    ;;
  h3n2)
    PARAMS=benchmarks/results/params.json
    TRAIN_DATA=data/h3n2/train
    CHECKPOINT="${CHECKPOINT:-checkpoints/h3n2_v3_lit_hotspot/best.pt}"
    MAX_SEQ_LEN=566
    PIN_GROUPS="${PIN_GROUPS:-2}"
    METHODS="${METHODS:-neutral_bd artreeformer_adapted treesbm}"
    ;;
  hiv_geo)
    PARAMS=benchmarks/results/params_hiv_geo.json
    TRAIN_DATA=data/hiv_geo/train
    CHECKPOINT="${CHECKPOINT:-checkpoints/hiv_geo_v1/best.pt}"
    MAX_SEQ_LEN=900
    PIN_GROUPS="${PIN_GROUPS:-}"
    # AR pools are flu/COVID-sized; Neutral+TreeSBM for HIV.
    METHODS="${METHODS:-neutral_bd treesbm}"
    ;;
  hiv_temporal)
    PARAMS=benchmarks/results/params_hiv_temporal.json
    TRAIN_DATA=data/hiv_temporal/train
    CHECKPOINT="${CHECKPOINT:-checkpoints/hiv_temporal_v1/best.pt}"
    MAX_SEQ_LEN=900
    PIN_GROUPS="${PIN_GROUPS:-}"
    METHODS="${METHODS:-neutral_bd treesbm}"
    ;;
  *)
    echo "Unknown VIRUS=$VIRUS"; exit 1
    ;;
esac

OUT="${OUT:-benchmarks/results/table_c3_${VIRUS}_N${N}.csv}"
PER_ROOT="${PER_ROOT:-benchmarks/results/table_c3_${VIRUS}_N${N}_per_root.csv}"
CUTS="${CUTS:-benchmarks/results/tables/table_c3_${VIRUS}_h_cuts.json}"
CACHE_DIR="${CACHE_DIR:-benchmarks/results/coverage_cache_c3_${VIRUS}_N${N}}"

echo "Start: $(date) virus=$VIRUS ckpt=$CHECKPOINT L=$MAX_SEQ_LEN pin=$PIN_GROUPS methods=$METHODS"
if [[ ! -f "$PARAMS" ]]; then
  $PYTHON benchmarks/fit_params.py --train-data "$TRAIN_DATA" --out "$PARAMS"
fi

PIN_ARGS=()
if [[ -n "${PIN_GROUPS}" ]]; then
  PIN_ARGS=(--pin-groups "$PIN_GROUPS")
fi

# shellcheck disable=SC2086
$PYTHON benchmarks/table_c3_horizon.py \
  --virus "$VIRUS" \
  --N "$N" \
  --K-list $K_LIST \
  --max-roots "$MAX_ROOTS" \
  --methods $METHODS \
  --checkpoint "$CHECKPOINT" \
  --max-seq-len "$MAX_SEQ_LEN" \
  --e-list "$E_LIST" \
  --n-steps 50 \
  --cache-dir "$CACHE_DIR" \
  --out "$OUT" \
  --per-root-out "$PER_ROOT" \
  --cuts-out "$CUTS" \
  "${PIN_ARGS[@]}"

echo "Done: $(date) -> $OUT cuts=$CUTS"
