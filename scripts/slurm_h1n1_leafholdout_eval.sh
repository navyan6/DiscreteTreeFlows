#!/bin/bash
#SBATCH --job-name=h1n1_lh_eval
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/h1n1_lh_eval_%j.log
#SBATCH --error=logs/h1n1_lh_eval_%j.log
#
# Standalone leaf-holdout recovery eval, split out from slurm_h1n1_leafholdout_train.sh
# so it can be re-run on its own against an already-trained checkpoint without
# retraining (the checkpoint from the first run is already at
# checkpoints/h1n1_leafholdout_v1/best.pt).

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}

cd ~/DiscreteTreeFlows
mkdir -p logs

CKPT=checkpoints/h1n1_leafholdout_v1/best.pt
if [[ ! -f "$CKPT" && -f "$LABHOME/checkpoints_backup/h1n1_leafholdout_v1/best.pt" ]]; then
    mkdir -p checkpoints/h1n1_leafholdout_v1
    cp -n "$LABHOME/checkpoints_backup/h1n1_leafholdout_v1/best.pt" "$CKPT"
fi

echo "Start: $(date)"

$PYTHON -u scripts/eval_leaf_holdout.py \
    --checkpoint "$CKPT" \
    --data       data/h1n1_leafholdout \
    --n-steps    30

echo "Done: $(date)"
