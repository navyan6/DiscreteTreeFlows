#!/bin/bash
#SBATCH --job-name=treesbm_train
#SBATCH --partition=gpu          # adjust to your Penn partition (gpu, a100, etc.)
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm_train_%j.log
#SBATCH --error=logs/slurm_train_%j.log

# ── environment ───────────────────────────────────────────────────────────────
# Replace <your-conda-env> with your environment name
conda activate <your-conda-env>

# ── paths ─────────────────────────────────────────────────────────────────────
# Run from repo root — adjust if submitting from elsewhere
cd /path/to/DiscreteTreeFlows        # ← set your project path on the cluster

mkdir -p logs checkpoints

echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURMD_NODENAME"
echo "GPU:    $CUDA_VISIBLE_DEVICES"
echo "Start:  $(date)"

# ── train ─────────────────────────────────────────────────────────────────────
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
