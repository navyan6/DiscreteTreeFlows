#!/bin/bash
#SBATCH --job-name=epidemic_pipe
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=48:00:00
#SBATCH --output=logs/epidemic_pipeline_%j.log
#SBATCH --error=logs/epidemic_pipeline_%j.log
#
# Generic tree pipeline for epidemic/outbreak split dirs.
# Env: DATA_ROOT, optional PREFIX_TRAIN / PREFIX_VAL / PREFIX_TEST

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

DATA_ROOT="${DATA_ROOT:?set DATA_ROOT}"
case "$DATA_ROOT" in
  *filo_l*) PT="${PREFIX_TRAIN:-filo_train}"; PV="${PREFIX_VAL:-filo_val}"; PX="${PREFIX_TEST:-filo_test}" ;;
  *covid_epidemic*) PT="${PREFIX_TRAIN:-covidetrain}"; PV="${PREFIX_VAL:-covideval}"; PX="${PREFIX_TEST:-covidetest}" ;;
  *h3n2_epidemic*) PT="${PREFIX_TRAIN:-h3n2etrain}"; PV="${PREFIX_VAL:-h3n2eval}"; PX="${PREFIX_TEST:-h3n2etest}" ;;
  *h1n1_epidemic*) PT="${PREFIX_TRAIN:-h1n1etrain}"; PV="${PREFIX_VAL:-h1n1eval}"; PX="${PREFIX_TEST:-h1n1etest}" ;;
  *panflu_forecast*|*h3n2_forecast*) PT="${PREFIX_TRAIN:-pflutrain}"; PV="${PREFIX_VAL:-pfluval}"; PX="${PREFIX_TEST:-pflutest}" ;;
  *) PT="${PREFIX_TRAIN:?set PREFIX_TRAIN}"; PV="${PREFIX_VAL:?set PREFIX_VAL}"; PX="${PREFIX_TEST:?set PREFIX_TEST}" ;;
esac

echo "Start: $(date) DATA_ROOT=$DATA_ROOT"

for pair in "train:${PT}" "val:${PV}" "test:${PX}"; do
  s="${pair%%:*}"
  p="${pair##*:}"
  dir="$DATA_ROOT/$s"
  if [[ ! -d "$dir" ]] || [[ -z "$(ls -A "$dir"/*_group_*.fasta 2>/dev/null || true)" ]]; then
    echo "SKIP $s (no groups in $dir)"
    continue
  fi
  echo "=== $s prefix=$p ==="
  $PYTHON scripts/run_all_groups.py --data-dir "$dir" --prefix "$p" \
    --workers 16 --stop-after refine
done

echo "Pipeline done: $(date)"
