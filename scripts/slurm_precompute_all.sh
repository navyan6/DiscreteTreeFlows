#!/bin/bash
#SBATCH --job-name=treesbm_precompute
#SBATCH --partition=dgx-b200
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=logs/slurm_precompute_%j.log
#SBATCH --error=logs/slurm_precompute_%j.log

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURMD_NODENAME"
echo "GPU:    $CUDA_VISIBLE_DEVICES"
echo "Start:  $(date)"

echo "=== PLM embeddings ==="
$PYTHON scripts/precompute_plm.py

echo "=== Reference rates ==="
$PYTHON scripts/precompute_ref_rates.py --max-seq-len 566 --batch-size 8

echo "Done: $(date)"
