#!/bin/bash
#SBATCH --job-name=covid_pipeline
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/covid_pipeline_%j.log
#SBATCH --error=logs/covid_pipeline_%j.log
#
# COVID Spike geographic-split pipeline (CPU): mafft -> fasttree -> augur refine
# (treetime root/timetree) -> augur ancestral (ASR) -> translate, per group.
#
# Prereq (run first):
#   sbatch scripts/slurm_covid_extract.sh              (nextclade Spike extraction)
#   python scripts/prepare_covid_geo.py                (country grouping + split)
#
# Resumable: run_all_groups skips completed stages, so re-submit if it times out.
# --stop-after translate skips the legacy GPU embed stage (train uses precompute_plm).

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Start: $(date)"

$PYTHON scripts/run_all_groups.py --data-dir data/covid/train --prefix covidtrain \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/covid/val   --prefix covidval \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/covid/test  --prefix covidtest \
    --workers 16 --stop-after translate

echo "Pipeline done: $(date)"
