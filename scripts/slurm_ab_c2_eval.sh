#!/bin/bash
#SBATCH --job-name=ab_c2_eval
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=4:00:00
#SBATCH --output=logs/ab_c2_eval_%j.log
#SBATCH --error=logs/ab_c2_eval_%j.log
#
# CPU eval: Cov@100/500 e2 + viral clade_recall / mean_min_edit on Rod.82 samples.
#
#   sbatch scripts/slurm_ab_c2_eval.sh
#   sbatch --dependency=afterok:A:B scripts/slurm_ab_c2_eval.sh

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p logs benchmarks/results/tables

OUT="${OUT:-benchmarks/results/tables/table_c2_ab_metrics.json}"
MODELS="${MODELS:-plm_prior,treesbm_ab_oas}"

echo "Start: $(date -Is) models=$MODELS out=$OUT"
$PYTHON -u scripts/eval_ab_c2_from_samples.py \
  --models "$MODELS" \
  --K-list 100,500 \
  --viral-K 100 \
  --out "$OUT" \
  --md-out "${OUT%.json}.md"

echo "Done: $(date -Is) -> $OUT"
