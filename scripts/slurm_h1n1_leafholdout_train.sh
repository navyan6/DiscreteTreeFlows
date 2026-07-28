#!/bin/bash
#SBATCH --job-name=h1n1_lh_train
#SBATCH --partition=b200-mig90
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/h1n1_lh_train_%j.log
#SBATCH --error=logs/h1n1_lh_train_%j.log
#
# H1N1 HA LEAF-HOLDOUT precompute + FRESH bridge-matched training (retrain from
# scratch -- own checkpoint dir checkpoints/h1n1_leafholdout_v1, no --resume, no
# reuse of any other checkpoint). Run AFTER slurm_h1n1_leafholdout_pipeline.sh.
#
# Split here is a RANDOM tree-level split (not geographic); the actual holdout
# under test is per-tree leaves, pulled out before each tree was even built by
# prepare_h1n1_leafholdout.py. After this finishes, run:
#   python scripts/eval_leaf_holdout.py --checkpoint checkpoints/h1n1_leafholdout_v1/best.pt \
#       --data data/h1n1_leafholdout
# to check whether the process recovers the real held-out leaves per tree.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints/h1n1_leafholdout_v1

echo "Start: $(date)"

# 1. ESM-2 embeddings + reference (ESM-only) mutation rates, per split
for split in train val test; do
    echo "=== precompute $split ==="
    $PYTHON scripts/precompute_plm.py       --data data/h1n1_leafholdout/$split
    $PYTHON scripts/precompute_ref_rates.py --data data/h1n1_leafholdout/$split --max-seq-len 566
done

# 2. Fresh bridge-matched training (no --resume: this is a from-scratch retrain)
$PYTHON -u scripts/train.py \
    --data        data/h1n1_leafholdout/train \
    --val-data    data/h1n1_leafholdout/val \
    --test-data   data/h1n1_leafholdout/test \
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
    --ckpt-dir    checkpoints/h1n1_leafholdout_v1

echo "Training done: $(date)"

# 3. Leaf-holdout recovery eval (generate from each tree's real root, compare
#    against the real leaves that were withheld before the tree was built)
$PYTHON -u scripts/eval_leaf_holdout.py \
    --checkpoint checkpoints/h1n1_leafholdout_v1/best.pt \
    --data       data/h1n1_leafholdout \
    --n-steps    30

echo "Done: $(date)"
