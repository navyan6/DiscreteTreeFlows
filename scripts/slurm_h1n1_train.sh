#!/bin/bash
#SBATCH --job-name=h1n1_train
#SBATCH --partition=b200-mig90
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/h1n1_train_%j.log
#SBATCH --error=logs/h1n1_train_%j.log
#
# H1N1 HA precompute + bridge-matched training WITH empirical column entropy
# (the "updated framework"). Run AFTER slurm_h1n1_pipeline.sh completes.
#
# H1N1 HA is ~566 aa like H3N2 (NOT spike's 1273), so --max-seq-len 566 here,
# not 1280. Own checkpoint dir (checkpoints/h1n1_v1); does not touch h3n2/covid.
#
# --entropy-source empirical: per-column Shannon entropy of the H1N1 train
# alignment computed once at startup (watch the "Computing empirical column
# entropy..." line), fed to the mutation head + L_mut reweight, saved in the
# checkpoint, and reused at generation.
#
# --resume is safe (own dir): no-op on the first run, resumes this entropy
# checkpoint on resubmit after the 16h wall. Never point --ckpt-dir elsewhere.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints/h1n1_v1

echo "Start: $(date)"

# 1. ESM-2 embeddings + reference (ESM-only) mutation rates, per split
for split in train val test; do
    echo "=== precompute $split ==="
    $PYTHON scripts/precompute_plm.py       --data data/h1n1/$split
    $PYTHON scripts/precompute_ref_rates.py --data data/h1n1/$split --max-seq-len 566
done

# 2. Bridge-matched training with empirical column entropy
$PYTHON -u scripts/train.py \
    --data        data/h1n1/train \
    --val-data    data/h1n1/val \
    --test-data   data/h1n1/test \
    --max-seq-len 566 \
    --epochs      100 \
    --patience    50 \
    --bridge-c    1.0 \
    --lambda-mut  5.0 \
    --use-site-entropy \
    --use-entropy-loss-weighting \
    --entropy-source empirical \
    --entropy-weight-alpha 3.0 \
    --entropy-weight-floor 1.0 \
    --ckpt-dir    checkpoints/h1n1_v1 \
    --resume

echo "Done: $(date)"
