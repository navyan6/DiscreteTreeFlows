#!/bin/bash
#SBATCH --job-name=treesbm_eval
#SBATCH --partition=b200-mig90
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=logs/slurm_eval_%j.log
#SBATCH --error=logs/slurm_eval_%j.log

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURMD_NODENAME"
echo "Start:  $(date)"

$PYTHON -u scripts/eval_single_tree.py \
    --checkpoint          checkpoints/best.pt \
    --data                data/train \
    --group               1 \
    --n-steps             50 \
    --max-leaves          400 \
    --max-seq-len         566 \
    --branch-rate-scale   6.0 \
    --mutation-rate-scale 0.04 \
    --pll-prune-threshold -3.0 \
    --seed                42

echo "Done: $(date)"
