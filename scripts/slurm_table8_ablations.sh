#!/bin/bash
#SBATCH --job-name=t8_abl
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/table8_ablate_%j.log
#SBATCH --error=logs/table8_ablate_%j.log
#
# Table 8 COVID ablation evals (gen-time flags on covid_v5_mutrec) + optional retrain.
#
# ABLATE=no_bridge|no_tree_ctx|no_seq_branch|no_bl_eval|no_internal_seqs|full|train_no_bl|no_lit_mask|no_entropy
# MODE=enrich (default mut recall + pLM NLL) | coverage (K=100 e-abs) | baselines (Tree-KL/Split-KL/W1)
#
#   ABLATE=no_entropy sbatch --qos=mig-max scripts/slurm_table8_ablations.sh
#   ABLATE=full MODE=coverage sbatch --qos=mig-max --time=24:00:00 scripts/slurm_table8_ablations.sh
#   ABLATE=no_bridge MODE=baselines sbatch --qos=mig-max --time=8:00:00 scripts/slurm_table8_ablations.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints benchmarks/results

ABLATE="${ABLATE:-no_bridge}"
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
BASE_MAX_ROOTS="${BASE_MAX_ROOTS:-20}"
BASE_K="${BASE_K:-20}"
BASE_M="${BASE_M:-20}"

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

if [[ "$ABLATE" == "train_no_bl" ]]; then
  echo "Train ablation: without branch-length head (lambda_br=0)"
  export LAMBDA_BR=0
  export CKPT_DIR=checkpoints/covid_ablate_no_bl
  # Re-exec train recipe in-process (this job IS the train job).
  exec env LAMBDA_BR=0 CKPT_DIR=checkpoints/covid_ablate_no_bl \
    bash scripts/slurm_covid_train_mut_recovery.sh
fi

CKPT_PATH=$(resolve_ckpt "$CKPT")
EXTRA=()
TAG="$ABLATE"
case "$ABLATE" in
  no_bridge) EXTRA+=(--ablate-bridge);;
  no_tree_ctx) EXTRA+=(--ablate-tree-context);;
  no_seq_branch) EXTRA+=(--branching-mode poisson_ref --ref-lambda 1.0);;
  no_bl_eval) EXTRA+=(--ablate-branch-length-head);;
  no_internal_seqs) EXTRA+=(--ablate-internal-node-seqs);;
  no_entropy) EXTRA+=(--ablate-site-entropy);;
  full) EXTRA+=();;
  no_fitness) EXTRA+=(--fitness-beta 0);;
  no_lit_mask)
    EXTRA+=(--no-lit-hotspot-mask)
    MASK=""
    ;;
  *) echo "Unknown ABLATE=$ABLATE"; exit 1;;
esac

echo "Start: $(date) ablate=$ABLATE mode=$MODE ckpt=$CKPT_PATH extra=${EXTRA[*]}"

if [[ "$MODE" == "enrich" ]]; then
  OUT="checkpoints/table8_${TAG}_mrs${MRS}.json"
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
  # Coverage@K e-abs (includes coverage_obs_e2 / e5). Focused treesbm-only.
  PARAMS=benchmarks/results/params_covid.json
  if [[ ! -f "$PARAMS" ]]; then
    $PYTHON benchmarks/fit_params.py --train-data "$TRAIN" --out "$PARAMS"
  fi
  OUT="benchmarks/results/coverage_curves_covid_N16_table8_${TAG}_K${K_MAX}.csv"
  CACHE_DIR="benchmarks/results/coverage_cache_covid_t8_${TAG}"
  # Drop --no-lit-hotspot-mask for coverage (generation-only flag).
  COV_EXTRA=()
  for a in "${EXTRA[@]}"; do
    [[ "$a" == "--no-lit-hotspot-mask" ]] && continue
    COV_EXTRA+=("$a")
  done
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
    "${COV_EXTRA[@]}"
  echo "Done: $(date) -> $OUT"
elif [[ "$MODE" == "baselines" ]]; then
  # Tree-KL / Split-KL / Branch W1 for this ablation only (treesbm, N=16).
  PARAMS=benchmarks/results/params_covid.json
  if [[ ! -f "$PARAMS" ]]; then
    $PYTHON benchmarks/fit_params.py --train-data "$TRAIN" --out "$PARAMS"
  fi
  OUT="benchmarks/results/results_baselines_covid_t8_${TAG}.csv"
  BASE_EXTRA=()
  for a in "${EXTRA[@]}"; do
    [[ "$a" == "--no-lit-hotspot-mask" ]] && continue
    BASE_EXTRA+=("$a")
  done
  $PYTHON benchmarks/run_table.py \
    --test-data "$DATA" \
    --train-data "$TRAIN" \
    --params "$PARAMS" \
    --N 16 --K "$BASE_K" --M "$BASE_M" --max-roots "$BASE_MAX_ROOTS" \
    --max-seq-len 1280 \
    --methods treesbm \
    --checkpoint "$CKPT_PATH" \
    --out "$OUT" \
    "${BASE_EXTRA[@]}"
  echo "Done: $(date) -> $OUT"
else
  echo "Unknown MODE=$MODE (enrich|coverage|baselines)"; exit 1
fi
