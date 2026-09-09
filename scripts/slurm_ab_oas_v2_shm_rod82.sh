#!/bin/bash
#SBATCH --job-name=ab_oas_v2_rod
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/ab_oas_v2_rod82_%j.log
#SBATCH --error=logs/ab_oas_v2_rod82_%j.log
#
# Rodriguez 82 eval for ab_oas_1m_v2_shm → samples/treesbm_ab_oas_v2
# (does not overwrite treesbm / treesbm_ab / treesbm_ab_oas).
#
#   sbatch --dependency=afterok:<train_job> scripts/slurm_ab_oas_v2_shm_rod82.sh

set -euo pipefail

export CKPT="${CKPT:-checkpoints/ab_oas_1m_v2_shm/best.pt}"
export WRITE_AS="${WRITE_AS:-treesbm_ab_oas_v2}"
export CONFIG="${CONFIG:-antibody_benchmark/configs/full_ab_oas_treesbm.yaml}"

# Guard
for keep in treesbm treesbm_ab treesbm_ab_oas; do
  if [[ "$WRITE_AS" == "$keep" ]]; then
    echo "REFUSING write-as=$WRITE_AS" >&2
    exit 2
  fi
done

exec bash scripts/slurm_ab_oas_rod82_eval.sh
