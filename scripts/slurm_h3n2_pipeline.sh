#!/bin/bash
#SBATCH --job-name=h3n2_pipeline
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/h3n2_pipeline_%j.log
#SBATCH --error=logs/h3n2_pipeline_%j.log
#
# H3N2 temporal pipeline (CPU): mafft -> fasttree -> augur refine (treetime root/
# timetree) -> augur ancestral (ASR) -> translate, for each temporal split into
# its OWN dir (no group_NNN collision with the existing multi-subtype data).
#
# Prereq (run locally, then rsync the split dirs up):
#   python scripts/prepare_h3n2_temporal.py
#   rsync -av data/h3n2/  <cluster>:~/DiscreteTreeFlows/data/h3n2/
#
# Resumable: run_all_groups skips completed stages, so re-submit if it times out.
# --stop-after translate skips the legacy GPU embed stage (train uses precompute_plm).

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Start: $(date)"

$PYTHON scripts/run_all_groups.py --data-dir data/h3n2/train --prefix h3n2train \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/h3n2/val   --prefix h3n2val \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/h3n2/test  --prefix h3n2test \
    --workers 16 --stop-after translate

echo "Pipeline done: $(date)"
