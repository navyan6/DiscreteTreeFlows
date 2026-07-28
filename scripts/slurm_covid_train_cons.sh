#!/bin/bash
#SBATCH --job-name=covid_train_cons
#SBATCH --partition=b200-mig90
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/covid_train_cons_%j.log
#SBATCH --error=logs/covid_train_cons_%j.log
#
# COVID Spike training with empirical entropy AND inverse-entropy L_cons
# weighting (--use-entropy-cons-weighting): cold/low-entropy conserved sites
# are penalized hardest for mutating, so the model learns to stop over-mutating
# conserved regions (retention was 0.65 at generation w/o it).
#
# Fresh run into checkpoints/covid_v3_cons -- a clean experiment alongside
# covid_v2_entropy (which has the entropy feature + L_mut weight but no L_cons
# weight). --resume is safe (own dir): no-op first run, resumes on timeout.
# Precompute is already done from covid_v2 (PLM + ref-rate caches shared across
# runs on the same data), so it's skipped here.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints/covid_v3_cons

echo "Start: $(date)"

$PYTHON -u scripts/train.py \
    --data        data/covid/train \
    --val-data    data/covid/val \
    --test-data   data/covid/test \
    --max-seq-len 1280 \
    --epochs      100 \
    --patience    50 \
    --bridge-c    1.0 \
    --lambda-mut  5.0 \
    --use-site-entropy \
    --use-entropy-loss-weighting \
    --use-entropy-cons-weighting \
    --entropy-source empirical \
    --entropy-weight-alpha 3.0 \
    --entropy-weight-floor 1.0 \
    --ckpt-dir    checkpoints/covid_v3_cons \
    --resume

echo "Done: $(date)"
