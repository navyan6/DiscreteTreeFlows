#!/bin/bash
#SBATCH --job-name=hiv_pipeline
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=48:00:00
#SBATCH --output=logs/hiv_pipeline_%j.log
#SBATCH --error=logs/hiv_pipeline_%j.log
#
# HIV-1 Env TreeSBM refine pipeline (CPU):
#   mafft -> fasttree -> augur refine (TreeTime numdate bl.json) ->
#   augur ancestral -> translate to AA
#
# Prereq:
#   rsync data/hiv/*.fasta to Betty
#   python scripts/prepare_hiv_splits.py
#
# Env:
#   SPLIT=temporal|geo|both   (default both)
#   WORKERS=16

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"

cd ~/DiscreteTreeFlows
mkdir -p logs

SPLIT="${SPLIT:-both}"
WORKERS="${WORKERS:-16}"

echo "Start: $(date -Is) host=$(hostname) job=${SLURM_JOB_ID:-local} SPLIT=$SPLIT"

run_one() {
  local data_dir="$1"
  local prefix="$2"
  if [[ ! -d "$data_dir" ]]; then
    echo "SKIP missing $data_dir"
    return 0
  fi
  local n
  n=$(ls "$data_dir"/${prefix}_group_*.fasta 2>/dev/null | wc -l | tr -d ' ')
  echo "=== $data_dir prefix=$prefix groups=$n ==="
  if [[ "$n" -eq 0 ]]; then
    echo "ERROR: no groups in $data_dir" >&2
    return 1
  fi
  $PYTHON scripts/run_all_groups.py \
    --data-dir "$data_dir" \
    --prefix "$prefix" \
    --workers "$WORKERS" \
    --stop-after translate
}

if [[ "$SPLIT" == "temporal" || "$SPLIT" == "both" ]]; then
  run_one data/hiv_temporal/train hivtemporaltrain
  run_one data/hiv_temporal/val   hivtemporalval
  run_one data/hiv_temporal/test  hivtemporaltest
fi

if [[ "$SPLIT" == "geo" || "$SPLIT" == "both" ]]; then
  run_one data/hiv_geo/train hivgeotrain
  run_one data/hiv_geo/val   hivgeoval
  run_one data/hiv_geo/test  hivgeotest
fi

echo "Pipeline done: $(date -Is)"
