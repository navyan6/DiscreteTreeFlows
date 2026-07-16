#!/bin/bash
#SBATCH --job-name=h3n2_train
#SBATCH --partition=b200-mig90
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/h3n2_train_%j.log
#SBATCH --error=logs/h3n2_train_%j.log
#
# H3N2 temporal precompute + bridge-matched training (GPU).
# Run AFTER slurm_h3n2_pipeline.sh completes (needs group_NNN_rooted.nwk /
# _anc_aa.fasta / _bl.json in each split dir).

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints/h3n2_temporal

echo "Start: $(date)"

# 1. ESM-2 embeddings + reference (ESM-only) mutation rates, per split
for split in train val test; do
    echo "=== precompute $split ==="
    $PYTHON scripts/precompute_plm.py       --data data/h3n2/$split
    $PYTHON scripts/precompute_ref_rates.py --data data/h3n2/$split
done

# 2. Bridge-matched training (new architecture) on the temporal split
$PYTHON -u scripts/train.py \
    --data       data/h3n2/train \
    --val-data   data/h3n2/val \
    --test-data  data/h3n2/test \
    --epochs     100 \
    --patience   50 \
    --bridge-c   1.0 \
    --lambda-mut 5.0 \
    --ckpt-dir   checkpoints/h3n2_temporal

echo "Done: $(date)"
