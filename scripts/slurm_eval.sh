#!/bin/bash
#SBATCH --job-name=treesbm_eval
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=logs/slurm_eval_%j.log
#SBATCH --error=logs/slurm_eval_%j.log

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURMD_NODENAME"
echo "Start:  $(date)"

$PYTHON -u scripts/eval_test_set.py \
    --checkpoint  checkpoints/best.pt \
    --data        data/train \
    --split       checkpoints/split_indices.json \
    --n-steps     50 \
    --max-trees   5 \
    --max-leaves  200 \
    --max-seq-len 566

echo "Done: $(date)"
