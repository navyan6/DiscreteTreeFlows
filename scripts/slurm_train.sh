#!/bin/bash
#SBATCH --job-name=treesbm_train
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=logs/slurm_train_%j.log
#SBATCH --error=logs/slurm_train_%j.log

source ~/.bashrc
module load miniconda3/25.5.1
source activate treesbm

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURMD_NODENAME"
echo "GPU:    $CUDA_VISIBLE_DEVICES"
echo "Start:  $(date)"

python scripts/train.py \
    --data        data/train \
    --epochs      300 \
    --lr          1e-4 \
    --t-max       0.95 \
    --lambda-top  0.1 \
    --lambda-br   0.1 \
    --ckpt-dir    checkpoints \
    --max-seq-len 566

echo "Done: $(date)"
