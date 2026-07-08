#!/bin/bash
#SBATCH --job-name=treesbm_eval
#SBATCH --partition=dgx-b200
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=logs/slurm_eval_%j.log
#SBATCH --error=logs/slurm_eval_%j.log
#SBATCH --exclude=dgx024

source ~/.bashrc
module load miniconda3/25.5.1
source activate treesbm

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURMD_NODENAME"
echo "GPU:    $CUDA_VISIBLE_DEVICES"
echo "Start:  $(date)"

python -u scripts/eval_test_set.py \
    --checkpoint checkpoints/best.pt \
    --data       data/train \
    --split      checkpoints/split_indices.json \
    --n-steps    30 \
    --max-seq-len 566

echo "Done: $(date)"
