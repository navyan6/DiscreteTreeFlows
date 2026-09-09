#!/bin/bash
#SBATCH --job-name=covid_clade_pipe
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/covid_cladeholdout_pipeline_%j.log
#SBATCH --error=logs/covid_cladeholdout_pipeline_%j.log
#
# COVID Spike CLADE-HOLDOUT pipeline (CPU): mafft -> fasttree -> augur refine
# -> augur ancestral -> translate on reduced trees (held-out clade leaves
# already removed by prepare_covid_cladeholdout.py).
#
# Prereq:
#   sbatch scripts/slurm_covid_extract.sh
#   python scripts/prepare_covid_cladeholdout.py --write-clade-tsv
#   # (or pass existing --clade-tsv; falls back to year-proxy without labels)
#
# Eval after train:
#   python scripts/eval_leaf_holdout.py \
#       --checkpoint checkpoints/covid_cladeholdout_v1/best.pt \
#       --data data/covid_cladeholdout --max-seq-len 1280
#
# Resumable: run_all_groups skips completed stages.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Start: $(date)"

$PYTHON scripts/run_all_groups.py --data-dir data/covid_cladeholdout/train --prefix covidchtrain \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/covid_cladeholdout/val   --prefix covidchval \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/covid_cladeholdout/test  --prefix covidchtest \
    --workers 16 --stop-after translate

echo "Pipeline done: $(date)"
