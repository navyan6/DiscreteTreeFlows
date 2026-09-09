#!/bin/bash
#SBATCH --job-name=ab_oas_v4_rod
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/ab_oas_v4_cosinehead_rod82_%j.log
#SBATCH --error=logs/ab_oas_v4_cosinehead_rod82_%j.log
#
# Rodriguez 82 eval for ab_oas_1m_v4_cosinehead → samples/treesbm_ab_oas_v4
# (does not overwrite treesbm / treesbm_ab / treesbm_ab_oas / _v2 / _v3).
#
#   sbatch --dependency=afterok:<train_job> scripts/slurm_ab_oas_v4_cosinehead_rod82.sh

set -euo pipefail

export CKPT="${CKPT:-checkpoints/ab_oas_1m_v4_cosinehead/best.pt}"
export WRITE_AS="${WRITE_AS:-treesbm_ab_oas_v4}"
export CONFIG="${CONFIG:-antibody_benchmark/configs/full_ab_oas_treesbm_v4.yaml}"
export OAS_RESULTS="${OAS_RESULTS:-antibody_benchmark/results_ab_oas_v4}"
export T6_MODELS="${T6_MODELS:-thrifty,dasm_thrifty,cosine,treesbm,treesbm_ab,treesbm_ab_oas,treesbm_ab_oas_v2,treesbm_ab_oas_v3,treesbm_ab_oas_v4}"

for keep in treesbm treesbm_ab treesbm_ab_oas treesbm_ab_oas_v2 treesbm_ab_oas_v3; do
  if [[ "$WRITE_AS" == "$keep" ]]; then
    echo "REFUSING write-as=$WRITE_AS" >&2
    exit 2
  fi
done

exec bash scripts/slurm_ab_oas_rod82_eval.sh
