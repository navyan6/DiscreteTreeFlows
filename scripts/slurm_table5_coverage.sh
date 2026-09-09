#!/bin/bash
#SBATCH --job-name=t5_cov
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/table5_coverage_%j.log
#SBATCH --error=logs/table5_coverage_%j.log
#
# Table 5 Coverage@100 + min dist for COVID / H1N1 / HIV Env (pLM / AR / TreeSBM).
#
# COVID:
#   VIRUS=covid sbatch --qos=mig-max scripts/slurm_table5_coverage.sh
# H1N1:
#   VIRUS=h1n1 sbatch --qos=mig-max scripts/slurm_table5_coverage.sh
# HIV temporal / geo (max_seq_len=900; TreeSBM-only by default):
#   VIRUS=hiv_temporal METHODS="treesbm" sbatch --qos=mig-max scripts/slurm_table5_coverage.sh
#   VIRUS=hiv_geo METHODS="treesbm" sbatch --qos=mig-max scripts/slurm_table5_coverage.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results

VIRUS="${VIRUS:-covid}"
K_MAX="${K_MAX:-100}"
K_STEP="${K_STEP:-10}"
MAX_ROOTS="${MAX_ROOTS:-5}"
N_LIST="${N_LIST:-16}"

if [[ "$VIRUS" == "covid" ]]; then
  TEST_DATA=data/covid/test
  TRAIN_DATA=data/covid/train
  PARAMS=benchmarks/results/params_covid.json
  # Primary Table 5 enrichment ckpt is covid_v5_mutrec (not v7). Morning job
  # 7539947 used v7_pmc_hotspot via the old default — override with CHECKPOINT=…
  CHECKPOINT="${CHECKPOINT:-checkpoints/covid_v5_mutrec/best.pt}"
  MAX_SEQ_LEN=1280
  OUT="${OUT:-benchmarks/results/coverage_curves_covid_N16_table5_eabs.csv}"
  CACHE_DIR="${CACHE_DIR:-benchmarks/results/coverage_cache_covid_N16_eabs}"
elif [[ "$VIRUS" == "h1n1" ]]; then
  TEST_DATA=data/h1n1/test
  TRAIN_DATA=data/h1n1/train
  PARAMS=benchmarks/results/params_h1n1.json
  CHECKPOINT="${CHECKPOINT:-checkpoints/h1n1_v2_lit_hotspot/best.pt}"
  MAX_SEQ_LEN=566
  OUT="${OUT:-benchmarks/results/coverage_curves_h1n1_N16_table5_eabs.csv}"
  CACHE_DIR="${CACHE_DIR:-benchmarks/results/coverage_cache_h1n1_N16}"
elif [[ "$VIRUS" == "hiv_temporal" ]]; then
  TEST_DATA=data/hiv_temporal/test
  TRAIN_DATA=data/hiv_temporal/train
  PARAMS=benchmarks/results/params_hiv_temporal.json
  CHECKPOINT="${CHECKPOINT:-checkpoints/hiv_temporal_v1/best.pt}"
  MAX_SEQ_LEN=900
  OUT="${OUT:-benchmarks/results/coverage_curves_hiv_temporal_N16_table5_eabs.csv}"
  CACHE_DIR="${CACHE_DIR:-benchmarks/results/coverage_cache_hiv_temporal_N16_eabs}"
elif [[ "$VIRUS" == "hiv_geo" ]]; then
  TEST_DATA=data/hiv_geo/test
  TRAIN_DATA=data/hiv_geo/train
  PARAMS=benchmarks/results/params_hiv_geo.json
  CHECKPOINT="${CHECKPOINT:-checkpoints/hiv_geo_v1/best.pt}"
  MAX_SEQ_LEN=900
  OUT="${OUT:-benchmarks/results/coverage_curves_hiv_geo_N16_table5_eabs.csv}"
  CACHE_DIR="${CACHE_DIR:-benchmarks/results/coverage_cache_hiv_geo_N16_eabs}"
elif [[ "$VIRUS" == "h3n2" ]]; then
  TEST_DATA=data/h3n2/test
  TRAIN_DATA=data/h3n2/train
  PARAMS=benchmarks/results/params_h3n2.json
  CHECKPOINT="${CHECKPOINT:-checkpoints/h3n2_v3_lit_hotspot/best.pt}"
  MAX_SEQ_LEN=566
  OUT="${OUT:-benchmarks/results/coverage_curves_h3n2_N16_table5_eabs.csv}"
  CACHE_DIR="${CACHE_DIR:-benchmarks/results/coverage_cache_h3n2_N16_eabs}"
else
  echo "Unknown VIRUS=$VIRUS"; exit 1
fi

# HIV: TreeSBM-only unless METHODS is exported (AR topology pools are flu/COVID).
if [[ "$VIRUS" == "hiv_temporal" || "$VIRUS" == "hiv_geo" ]]; then
  METHODS="${METHODS:-treesbm}"
else
  METHODS="${METHODS:-neutral_bd plm_prior artreeformer_adapted treesbm}"
fi

echo "Start: $(date) virus=$VIRUS ckpt=$CHECKPOINT max_seq_len=$MAX_SEQ_LEN out=$OUT"
if [[ ! -f "$PARAMS" ]]; then
  $PYTHON benchmarks/fit_params.py --train-data "$TRAIN_DATA" --out "$PARAMS"
fi

# Paper Coverage for Table 5 should be absolute Hamming coverage_obs_e2
# (same e-abs protocol as H3N2 job 7459412). eps_frac=0.02 is legacy only.
# Prefer CHECKPOINT override: COVID best enrichment is covid_v5_mutrec unless
# v7 is explicitly selected.
E_LIST="${E_LIST:-0,1,2,3,5,8,10}"
EPS_FRAC="${EPS_FRAC:-0.02}"
SITE_TEMPERATURE="${SITE_TEMPERATURE:-1.0}"
EXTRA_ARGS=()
if [[ -n "${MAX_LEAVES:-}" ]]; then
  EXTRA_ARGS+=(--max-leaves "$MAX_LEAVES")
fi

# shellcheck disable=SC2086
$PYTHON benchmarks/coverage_curves.py \
  --test-data "$TEST_DATA" \
  --train-data "$TRAIN_DATA" \
  --params "$PARAMS" \
  --N $N_LIST \
  --K-max "$K_MAX" --K-step "$K_STEP" \
  --max-roots "$MAX_ROOTS" \
  --methods $METHODS \
  --eps-frac "$EPS_FRAC" \
  --e-list "$E_LIST" \
  --checkpoint "$CHECKPOINT" \
  --n-steps 50 \
  --max-seq-len "$MAX_SEQ_LEN" \
  --site-temperature "$SITE_TEMPERATURE" \
  --cache-dir "$CACHE_DIR" \
  --out "$OUT" \
  "${EXTRA_ARGS[@]}"

echo "Done: $(date) -> $OUT"
