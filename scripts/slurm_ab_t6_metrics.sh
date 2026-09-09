#!/bin/bash
#SBATCH --job-name=ab_t6_metrics
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=02:00:00
#SBATCH --output=logs/ab_t6_metrics_%j.log
#SBATCH --error=logs/ab_t6_metrics_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$HOME/DiscreteTreeFlows}"
mkdir -p logs benchmarks/results/tables

PY="${PY:-/vast/home/n/nnori/.conda/envs/treesbm/bin/python}"
[ -x "$PY" ] || PY="${PY_FALLBACK:-/vast/projects/pranam/lab/nnori/.conda/envs/treesbm/bin/python}"

echo "host=$(hostname) job=${SLURM_JOB_ID:-local} date=$(date -Is) py=$PY"
"$PY" -u scripts/eval_ab_t6_from_samples.py \
  --trees antibody_benchmark/data/processed/benchmark_trees.jsonl \
  --samples-dir antibody_benchmark/results/samples \
  --models thrifty,dasm_thrifty,cosine,treesbm,treesbm_ab,treesbm_ab_oas \
  --K 100 \
  --eps 1,2,3,5 \
  --out benchmarks/results/tables/table6_ab_track_c.json \
  --md-out benchmarks/results/tables/table6_ab_track_c.md

echo "DONE $(date -Is)"
