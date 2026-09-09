#!/bin/bash
#SBATCH --job-name=covid_temporal_pipe
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/covid_temporal_pipeline_%j.log
#SBATCH --error=logs/covid_temporal_pipeline_%j.log
#
# COVID Spike TEMPORAL-split pipeline (CPU): mafft -> fasttree -> augur refine
# -> augur ancestral -> translate under data/covid_temporal/.
#
# Prereq:
#   sbatch scripts/slurm_covid_extract.sh
#   python scripts/prepare_covid_temporal.py
#   # defaults: train≤2022 / val=2023 / test=2024–2025
#
# Resumable: run_all_groups skips completed stages.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Start: $(date)"

$PYTHON scripts/run_all_groups.py --data-dir data/covid_temporal/train --prefix covidttrain \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/covid_temporal/val   --prefix covidtval \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/covid_temporal/test  --prefix covidttest \
    --workers 16 --stop-after translate

echo "Pipeline done: $(date)"
