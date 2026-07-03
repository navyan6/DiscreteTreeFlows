#!/bin/bash
#SBATCH --job-name=treesbm_plm
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm_plm_%j.log
#SBATCH --error=logs/slurm_plm_%j.log

# Only needed if the data/train/*_plm.pt files aren't synced from local.
# The local precompute already ran on CPU — syncing those files is faster.
# Use this script only to regenerate caches on the cluster from scratch.

conda activate <your-conda-env>
cd /path/to/DiscreteTreeFlows

mkdir -p logs

echo "Job ID: $SLURM_JOB_ID  GPU: $CUDA_VISIBLE_DEVICES  Start: $(date)"
python scripts/precompute_plm.py
echo "Done: $(date)"
