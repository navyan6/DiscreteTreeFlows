#!/bin/bash
#SBATCH --job-name=h3n2_forecast_pipe
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/h3n2_temporal_forecast_pipeline_%j.log
#SBATCH --error=logs/h3n2_temporal_forecast_pipeline_%j.log
#
# H3N2 temporal-forecast pipeline for data/h3n2_temporal_forecast/ (alternate
# year cutoffs). For the default data/h3n2 protocol use slurm_h3n2_pipeline.sh.
#
# Prereq:
#   python scripts/prepare_h3n2_temporal.py \
#       --out-base data/h3n2_temporal_forecast \
#       --train-end-year 2023 --val-year 2024 \
#       --test-start-year 2025 --test-end-year 2025
#   # Note: test 2025 needs a 2025 source window; currently source windows
#   # top out at 2024 — extend SOURCE_WINDOWS when new FASTAs arrive.
#
# Resumable: run_all_groups skips completed stages.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

OUT=${1:-data/h3n2_temporal_forecast}

echo "Start: $(date)  OUT=$OUT"

$PYTHON scripts/run_all_groups.py --data-dir "$OUT/train" --prefix h3n2train \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir "$OUT/val"   --prefix h3n2val \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir "$OUT/test"  --prefix h3n2test \
    --workers 16 --stop-after translate

echo "Pipeline done: $(date)"
