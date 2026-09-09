#!/bin/bash
#SBATCH --job-name=h1n1_temporal_pipe
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/h1n1_temporal_pipeline_%j.log
#SBATCH --error=logs/h1n1_temporal_pipeline_%j.log
#
# H1N1 HA TEMPORAL-split pipeline (CPU): mafft -> fasttree -> augur refine
# -> augur ancestral -> translate, per group under data/h1n1_temporal/.
#
# Prereq:
#   python scripts/prepare_h1n1_temporal.py
#   # defaults: train≤2023 / val=2024 / test=2025
#   rsync -avP data/h1n1_temporal <cluster>:~/DiscreteTreeFlows/data/
#
# Resumable: run_all_groups skips completed stages.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Start: $(date)"

$PYTHON scripts/run_all_groups.py --data-dir data/h1n1_temporal/train --prefix h1n1ttrain \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/h1n1_temporal/val   --prefix h1n1tval \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/h1n1_temporal/test  --prefix h1n1ttest \
    --workers 16 --stop-after translate

echo "Pipeline done: $(date)"
