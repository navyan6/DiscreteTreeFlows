#!/bin/bash
#SBATCH --job-name=h1n1_lh_pipeline
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/h1n1_lh_pipeline_%j.log
#SBATCH --error=logs/h1n1_lh_pipeline_%j.log
#
# H1N1 HA LEAF-HOLDOUT split pipeline (CPU): mafft -> fasttree -> augur refine
# -> augur ancestral -> translate, per group. Same stages as slurm_h1n1_pipeline.sh
# but on data/h1n1_leafholdout/{train,val,test} (random tree-level split, each
# tree already missing its held-out leaves -- see prepare_h1n1_leafholdout.py).
#
# Prereq (run first):
#   python scripts/prepare_h1n1_leafholdout.py     (random split + per-tree leaf holdout)
#   rsync -avP data/h1n1_leafholdout <cluster>:~/DiscreteTreeFlows/data/
#
# Resumable: run_all_groups skips completed stages, so re-submit if it times out.
# --stop-after translate skips the legacy GPU embed stage (train uses precompute_plm).

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Start: $(date)"

$PYTHON scripts/run_all_groups.py --data-dir data/h1n1_leafholdout/train --prefix h1n1lhtrain \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/h1n1_leafholdout/val   --prefix h1n1lhval \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/h1n1_leafholdout/test  --prefix h1n1lhtest \
    --workers 16 --stop-after translate

echo "Pipeline done: $(date)"
