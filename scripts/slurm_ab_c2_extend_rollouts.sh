#!/bin/bash
#SBATCH --job-name=ab_c2_roll
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --time=12:00:00
#SBATCH --output=logs/ab_c2_rollouts_%j.log
#SBATCH --error=logs/ab_c2_rollouts_%j.log
#
# Extend Rod.82 rollouts so gen-leaf pool can support Cov@500 (median ~10 leaves/rollout
# → need ~50–100 rollouts). Skips existing rollout_*.json.
#
#   MODEL=plm_prior N_ROLLOUTS=100 sbatch --qos=mig-max scripts/slurm_ab_c2_extend_rollouts.sh
#   MODEL=treesbm WRITE_AS=treesbm_ab_oas CONFIG=antibody_benchmark/configs/full_ab_oas_treesbm.yaml \
#     N_ROLLOUTS=100 sbatch --qos=mig-max scripts/slurm_ab_c2_extend_rollouts.sh

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p logs antibody_benchmark/results/samples

MODEL="${MODEL:-plm_prior}"
WRITE_AS="${WRITE_AS:-}"
N_ROLLOUTS="${N_ROLLOUTS:-100}"
SAMPLES_DIR="${SAMPLES_DIR:-antibody_benchmark/results/samples}"

if [[ "$MODEL" == "treesbm" || "$MODEL" == "treesbm_ab_oas" ]]; then
  CONFIG="${CONFIG:-antibody_benchmark/configs/full_ab_oas_treesbm.yaml}"
  ROLL_MODEL=treesbm
  WRITE_AS="${WRITE_AS:-treesbm_ab_oas}"
else
  CONFIG="${CONFIG:-antibody_benchmark/configs/full.yaml}"
  ROLL_MODEL=plm_prior
  WRITE_AS=""
fi

echo "Start: $(date -Is) model=$ROLL_MODEL write_as=${WRITE_AS:-none} n_rollouts=$N_ROLLOUTS"
EXTRA=()
if [[ -n "$WRITE_AS" ]]; then
  EXTRA+=(--write-as "$WRITE_AS")
fi

$PYTHON antibody_benchmark/scripts/run_rollouts.py \
  --config "$CONFIG" \
  --models "$ROLL_MODEL" \
  --n-rollouts "$N_ROLLOUTS" \
  --out-dir "$SAMPLES_DIR" \
  "${EXTRA[@]}"

TARGET="${WRITE_AS:-$ROLL_MODEL}"
nsamp=$(find "$SAMPLES_DIR/$TARGET" -name 'rollout_*.json' 2>/dev/null | wc -l | tr -d ' ')
echo "samples_written_or_kept=$nsamp under $SAMPLES_DIR/$TARGET"
echo "Done: $(date -Is)"
