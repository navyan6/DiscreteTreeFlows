#!/bin/bash
#SBATCH --job-name=precompute_ds
#SBATCH --partition=b200-mig90
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=8G
#SBATCH --time=8:00:00
#SBATCH --output=logs/precompute_%x_%j.log
#SBATCH --error=logs/precompute_%x_%j.log
#
# Usage:
#   sbatch --qos=mig-max --job-name=pre_h1n1t \
#     --export=ALL,DATA_ROOT=data/h1n1_temporal,MAX_SEQ_LEN=566 \
#     scripts/_ops/slurm_precompute_dataset.sh

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
cd ~/DiscreteTreeFlows
mkdir -p logs

DATA_ROOT="${DATA_ROOT:?set DATA_ROOT}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"
BATCH_SIZE="${BATCH_SIZE:-8}"

echo "Job $SLURM_JOB_ID node=$SLURMD_NODENAME gpu=$CUDA_VISIBLE_DEVICES"
echo "DATA_ROOT=$DATA_ROOT MAX_SEQ_LEN=$MAX_SEQ_LEN Start: $(date)"

for split in train val test; do
  d="$DATA_ROOT/$split"
  echo "=== PLM $d ==="
  $PYTHON scripts/precompute_plm.py --data "$d"
  echo "=== REF_RATES $d ==="
  $PYTHON scripts/precompute_ref_rates.py --data "$d" --max-seq-len "$MAX_SEQ_LEN" --batch-size "$BATCH_SIZE"
done

echo "Done: $(date)"
