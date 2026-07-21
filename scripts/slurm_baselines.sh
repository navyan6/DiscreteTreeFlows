#!/bin/bash
#SBATCH --job-name=baselines
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --output=logs/baselines_%j.log
#SBATCH --error=logs/baselines_%j.log
#
# Baseline-only benchmark grid (neutral_bd, empirical_bd, plm_prior --
# TreeSBM excluded, ~170s/sample makes the full grid infeasible for now).

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results

echo "Start: $(date)"

$PYTHON benchmarks/run_table.py \
    --test-data data/h3n2/test \
    --params    benchmarks/results/params.json \
    --N 16 32 64 --K 100 --M 50 --max-roots 100 \
    --out benchmarks/results/results_baselines.csv

echo "Done: $(date)"
