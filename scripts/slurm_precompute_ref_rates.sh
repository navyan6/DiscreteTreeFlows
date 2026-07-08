#!/bin/bash
#SBATCH --job-name=treesbm_ref_rates
#SBATCH --partition=dgx-b200
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=logs/slurm_ref_rates_%j.log
#SBATCH --error=logs/slurm_ref_rates_%j.log

source ~/.bashrc
module load miniconda3/25.5.1
source activate treesbm

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURMD_NODENAME"
echo "GPU:    $CUDA_VISIBLE_DEVICES"
echo "Start:  $(date)"

python scripts/precompute_ref_rates.py \
    --data        data/train \
    --max-seq-len 566 \
    --batch-size  8

echo "Done: $(date)"
