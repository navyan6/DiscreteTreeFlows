#!/bin/bash
#SBATCH --job-name=covid_train_ent
#SBATCH --partition=b200-mig90
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/covid_train_entropy_%j.log
#SBATCH --error=logs/covid_train_entropy_%j.log
#
# COVID Spike bridge-matched training WITH empirical column entropy, into a
# SEPARATE checkpoint dir (checkpoints/covid_v2_entropy). Runs alongside
# checkpoints/covid_v1 (the non-entropy baseline) rather than clobbering it,
# so the two can be compared.
#
# --resume is SAFE here even though --use-site-entropy adds a site_entropy_proj
# param: this script only ever resumes its OWN dir (checkpoints/covid_v2_entropy),
# which is always an entropy checkpoint with matching architecture + optimizer
# state. On the first run best.pt doesn't exist, so train.py's --resume is a
# no-op and it starts fresh (site_entropy_proj is zero-initialized -> epoch 1
# begins identical to a non-entropy model, then learns the entropy correction).
# On a resubmit after the 16h wall (likely, at ~21 min/epoch) it picks up where
# it left off instead of restarting. Do NOT point --ckpt-dir at covid_v1 -- that
# non-entropy checkpoint's optimizer state would not match this architecture.
# col_entropy is recomputed deterministically at every startup, so resume stays
# consistent.
#
# --entropy-source empirical: Shannon entropy of each TRAIN-alignment column is
# computed once at startup (a ~1-2 min pass over all train groups; watch for the
# "Computing empirical column entropy..." line + its mean/max/nonzero summary),
# saved in the checkpoint, and reused at generation so train/inference match.
#
# --entropy-weight-alpha/floor 1.0 -> hotspot sites get up to 2x weight in L_mut
# (weight = floor + alpha*entropy, entropy in [0,1]). Raise alpha to 2-3 to make
# the model focus harder on antigenic (high-entropy) sites.
#
# Precompute (PLM + ref rates) is assumed already done by slurm_covid_train.sh;
# if not, uncomment the precompute loop below.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints/covid_v2_entropy

echo "Start: $(date)"

# # 1. (only if not already precomputed) ESM-2 embeddings + reference rates
# for split in train val test; do
#     echo "=== precompute $split ==="
#     $PYTHON scripts/precompute_plm.py       --data data/covid/$split
#     $PYTHON scripts/precompute_ref_rates.py --data data/covid/$split --max-seq-len 1280
# done

# 2. Bridge-matched training WITH empirical column entropy.
# --resume is a no-op on the first run (no best.pt yet) and picks up an entropy
# checkpoint on resubmit -- see header. Never repoints at covid_v1.
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
    --entropy-source empirical \
    --entropy-weight-alpha 3.0 \
    --entropy-weight-floor 1.0 \
    --ckpt-dir    checkpoints/covid_v2_entropy \
    --resume

echo "Done: $(date)"
