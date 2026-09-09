#!/bin/bash
#SBATCH --job-name=bdbv_pipeline
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/bdbv_pipeline_%j.log
#SBATCH --error=logs/bdbv_pipeline_%j.log
#
# Prereq:
#   python scripts/prepare_bdbv_temporal.py [--pan-ebolavirus]
#   rsync data/bdbv_temporal/ to cluster
#
# Env:
#   DATA_ROOT=data/bdbv_temporal  (default)
#   PREFIX_TRAIN=bdbvttrain PREFIX_VAL=bdbvtval PREFIX_TEST=bdbvttest

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

DATA_ROOT="${DATA_ROOT:-data/bdbv_temporal}"
if [[ "$DATA_ROOT" == *filo_l* ]]; then
  PT="${PREFIX_TRAIN:-filo_train}"
  PV="${PREFIX_VAL:-filo_val}"
  PX="${PREFIX_TEST:-filo_test}"
else
  PT="${PREFIX_TRAIN:-bdbvttrain}"
  PV="${PREFIX_VAL:-bdbvtval}"
  PX="${PREFIX_TEST:-bdbvttest}"
fi

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
    --workers 16 --stop-after translate
done

echo "Pipeline done: $(date)"
