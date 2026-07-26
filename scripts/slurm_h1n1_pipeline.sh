#!/bin/bash
#SBATCH --job-name=h1n1_pipeline
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/h1n1_pipeline_%j.log
#SBATCH --error=logs/h1n1_pipeline_%j.log
#
# H1N1 HA geolocation-split pipeline (CPU): mafft -> fasttree -> augur refine
# (treetime root/timetree) -> augur ancestral (ASR) -> translate, per group.
# Same pipeline as COVID; H1N1 HA is already a nucleotide CDS so there is NO
# gene-extraction step (unlike COVID spike).
#
# Prereq (run first):
#   # get the raw regional dumps onto the cluster (data/ is gitignored):
#   rsync -avP data/h1n1/train/*_h1n1.fasta <cluster>:~/DiscreteTreeFlows/data/h1n1/train/
#   python scripts/prepare_h1n1_geo.py     (geolocation grouping + 80/10/10 geo split)
#
# Resumable: run_all_groups skips completed stages, so re-submit if it times out.
# --stop-after translate skips the legacy GPU embed stage (train uses precompute_plm).

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Start: $(date)"

$PYTHON scripts/run_all_groups.py --data-dir data/h1n1/train --prefix h1n1train \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/h1n1/val   --prefix h1n1val \
    --workers 16 --stop-after translate
$PYTHON scripts/run_all_groups.py --data-dir data/h1n1/test  --prefix h1n1test \
    --workers 16 --stop-after translate

echo "Pipeline done: $(date)"
