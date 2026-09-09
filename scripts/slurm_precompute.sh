#!/bin/bash
#SBATCH --job-name=treesbm_plm
#SBATCH --partition=dgx-b200
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=logs/slurm_plm_%j.log
#SBATCH --error=logs/slurm_plm_%j.log


source ~/.bashrc
module load miniconda3/25.5.1
source activate treesbm

cd ~/DiscreteTreeFlows

mkdir -p logs

echo "Job ID: $SLURM_JOB_ID  GPU: $CUDA_VISIBLE_DEVICES  Start: $(date)"
python scripts/precompute_plm.py
echo "Done: $(date)"
