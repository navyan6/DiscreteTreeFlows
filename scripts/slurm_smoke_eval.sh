#!/bin/bash
#SBATCH --job-name=smoke_eval
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/smoke_eval_%j.log
#SBATCH --error=logs/smoke_eval_%j.log
#
# Quick single-tree generation + eval smoke test on a trained checkpoint (GPU).
# One tree, so a small mig45 slice + 1h is plenty; loads the saved col_entropy
# from the checkpoint automatically (eval_single_tree fix), so entropy models
# generate with the same signal they trained on.
#
# Usage: sbatch scripts/slurm_smoke_eval.sh <checkpoint> <data_dir> <max_seq_len> [group] [mut_rate_scale] [n_steps]
#   covid:              sbatch scripts/slurm_smoke_eval.sh checkpoints/covid_v2_entropy/best.pt data/covid/test 1280 1
#   sweep mut rate 0.3: sbatch scripts/slurm_smoke_eval.sh checkpoints/covid_v2_entropy/best.pt data/covid/test 1280 1 0.3
#   +deeper (200 steps): sbatch scripts/slurm_smoke_eval.sh checkpoints/covid_v2_entropy/best.pt data/covid/test 1280 1 0.3 200

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

CKPT="${1:?usage: sbatch slurm_smoke_eval.sh <checkpoint> <data_dir> <max_seq_len> [group] [mut_rate_scale] [n_steps]}"
DATA="${2:?usage: sbatch slurm_smoke_eval.sh <checkpoint> <data_dir> <max_seq_len> [group] [mut_rate_scale] [n_steps]}"
MSL="${3:?usage: sbatch slurm_smoke_eval.sh <checkpoint> <data_dir> <max_seq_len> [group] [mut_rate_scale] [n_steps]}"
GROUP="${4:-1}"
MRS="${5:-1.0}"      # --mutation-rate-scale (lower => fewer mutations; try 0.3 / 0.1 for over-mutation)
NSTEPS="${6:-50}"    # --n-steps (more => deeper trees)

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Start: $(date)  ckpt=$CKPT  data=$DATA  max_seq_len=$MSL  group=$GROUP  mut_rate_scale=$MRS  n_steps=$NSTEPS"
$PYTHON scripts/eval_single_tree.py \
    --checkpoint "$CKPT" \
    --data "$DATA" \
    --group "$GROUP" \
    --max-seq-len "$MSL" \
    --max-leaves 300 \
    --mutation-rate-scale "$MRS" \
    --n-steps "$NSTEPS"
echo "Done: $(date)"
