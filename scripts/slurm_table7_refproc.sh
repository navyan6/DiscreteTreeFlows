#!/bin/bash
#SBATCH --job-name=t7_ref
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/table7_refproc_%j.log
#SBATCH --error=logs/table7_refproc_%j.log
#
# Table 7 — Reference process ablation on held-out SARS-CoV-2 Brazil geo-test
# (data/covid/test). Gen-time R0 swap on covid_v5_mutrec.
#
# Rows (paper):
#   jtt_nofit         JTT/WAG/LG substitution, fitness no  (uses jtt)
#   esm2_650m_nofit   ESM-2-650M, fitness no
#   esm2_650m_fit     ESM-2-650M, fitness yes
#   esmc_nofit        ESM-C, fitness no
#   esmc_fit          ESM-C, fitness yes
#   progen2_fit       ProGen2-small (enijkamp), fitness yes
#   treesbm_fit       TreeSBM default (ESM-2-8M), fitness yes
#
# Rows 1–6 use --ablate-bridge (pure R0; no learned c_θ).
# Row 7 is full TreeSBM + --fitness-beta 1.0.
#
# MODE=enrich (default, mut recall) | coverage (K_MAX, default 100; set 1000 if feasible)
#
#   ROW=jtt_nofit sbatch --qos=mig-max --job-name=t7_jtt scripts/slurm_table7_refproc.sh
#   ROW=treesbm_fit MODE=coverage K_MAX=100 sbatch --qos=mig-max --time=24:00:00 \
#     --job-name=t7_cov_tsbm scripts/slurm_table7_refproc.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"
export PROGEN2_HOME="${PROGEN2_HOME:-$LABHOME/progen2}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints benchmarks/results

ROW="${ROW:-jtt_nofit}"
MODE="${MODE:-enrich}"
CKPT="${CKPT:-checkpoints/covid_v5_mutrec/best.pt}"
DATA="${DATA:-data/covid/test}"
TRAIN="${TRAIN:-data/covid/train}"
EVESCAPE=data/covid/evescape_spike_rbd.pt
MASK=results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt
MRS="${MRS:-0.5}"
NSTEPS="${NSTEPS:-100}"
MAX_TREES="${MAX_TREES:-20}"
K_MAX="${K_MAX:-100}"
K_STEP="${K_STEP:-10}"
MAX_ROOTS="${MAX_ROOTS:-5}"

resolve_ckpt() {
  local c="$1"
  if [[ -f "$c" ]]; then echo "$c"; return; fi
  local name; name=$(basename "$(dirname "$c")")
  local backup="$LABHOME/checkpoints_backup/${name}/best.pt"
  if [[ -f "$backup" ]]; then
    mkdir -p "$(dirname "$c")"
    cp -n "$backup" "$c" 2>/dev/null || cp "$backup" "$c"
    echo "$c"; return
  fi
  echo "$c"
}

CKPT_PATH=$(resolve_ckpt "$CKPT")
EXTRA=()
case "$ROW" in
  jtt_nofit|wag_nofit|lg_nofit)
    backend="${ROW%%_*}"; EXTRA+=(--r0-backend "$backend" --fitness-beta 0 --ablate-bridge);;
  esm2_650m_nofit) EXTRA+=(--r0-backend esm2_650m --fitness-beta 0 --ablate-bridge);;
  esm2_650m_fit)   EXTRA+=(--r0-backend esm2_650m --fitness-beta 1.0 --ablate-bridge);;
  esmc_nofit)      EXTRA+=(--r0-backend esmc --fitness-beta 0 --ablate-bridge);;
  esmc_fit)        EXTRA+=(--r0-backend esmc --fitness-beta 1.0 --ablate-bridge);;
  progen2_fit)     EXTRA+=(--r0-backend progen2 --fitness-beta 1.0 --ablate-bridge);;
  treesbm_fit)     EXTRA+=(--fitness-beta 1.0);;
  *) echo "Unknown ROW=$ROW"; exit 1;;
esac

# ESM-C: HF AutoModelForMaskedLM (Synthyra/ESMplusplus_small / biohub/ESMC-300M).
# Do NOT set PYTHONPATH=python_esmc — that shadows transformers → flash_attn/libcudart.
# ProGen2: official enijkamp/progen2 under $PROGEN2_HOME (progen2-small checkpoint).
if [[ "$ROW" == esmc_* ]]; then
  echo "ESM-C via HF (no PYTHONPATH side-install); PROGEN2_HOME unused"
fi
if [[ "$ROW" == progen2_* ]]; then
  echo "PROGEN2_HOME=$PROGEN2_HOME"
  if [[ ! -d "$PROGEN2_HOME/checkpoints/progen2-small" ]]; then
    echo "WARN: missing $PROGEN2_HOME/checkpoints/progen2-small"
  fi
fi

echo "Start: $(date) row=$ROW mode=$MODE ckpt=$CKPT_PATH extra=${EXTRA[*]}"

if [[ "$MODE" == "enrich" ]]; then
  OUT="checkpoints/table7_${ROW}_mrs${MRS}.json"
  $PYTHON scripts/eval_evescape_enrichment.py \
    --checkpoint "$CKPT_PATH" \
    --data "$DATA" \
    --max-seq-len 1280 \
    --evescape "$EVESCAPE" \
    --lit-hotspot-mask "$MASK" \
    --mutation-rate-scale "$MRS" \
    --n-steps "$NSTEPS" \
    --max-trees "$MAX_TREES" \
    --out "$OUT" \
    "${EXTRA[@]}"
  echo "Done: $(date) -> $OUT"
elif [[ "$MODE" == "coverage" ]]; then
  PARAMS=benchmarks/results/params_covid.json
  if [[ ! -f "$PARAMS" ]]; then
    $PYTHON benchmarks/fit_params.py --train-data "$TRAIN" --out "$PARAMS"
  fi
  OUT="benchmarks/results/coverage_curves_covid_N16_table7_${ROW}_K${K_MAX}.csv"
  CACHE_DIR="benchmarks/results/coverage_cache_covid_t7_${ROW}"
  # K=1000 is ~10× K=100; default K_MAX=100 unless explicitly overridden.
  $PYTHON benchmarks/coverage_curves.py \
    --test-data "$DATA" \
    --train-data "$TRAIN" \
    --params "$PARAMS" \
    --N 16 \
    --K-max "$K_MAX" --K-step "$K_STEP" \
    --max-roots "$MAX_ROOTS" \
    --methods treesbm \
    --eps-frac 0.02 \
    --e-list 0,1,2,3,5,8,10 \
    --checkpoint "$CKPT_PATH" \
    --n-steps 50 \
    --max-seq-len 1280 \
    --cache-dir "$CACHE_DIR" \
    --out "$OUT" \
    "${EXTRA[@]}"
  echo "Done: $(date) -> $OUT"
else
  echo "Unknown MODE=$MODE (enrich|coverage)"; exit 1
fi
