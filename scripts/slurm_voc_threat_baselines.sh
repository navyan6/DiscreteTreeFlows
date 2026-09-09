#!/bin/bash
#SBATCH --job-name=voc_bl
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/voc_bl_%j.log
#SBATCH --error=logs/voc_bl_%j.log
#
# NeutralBD (+ optional plm_prior) baselines for one VOC case.
#
#   CASE_DIR=results/voc_threat_panel/cases/Gamma_g003_test VOC_ID=Gamma \
#     sbatch scripts/slurm_voc_threat_baselines.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}

cd ~/DiscreteTreeFlows
mkdir -p logs

CASE_DIR="${CASE_DIR:?}"
VOC_ID="${VOC_ID:?}"
N_LEAVES="${N_LEAVES:-250}"
METHODS="${METHODS:-neutral_bd}"
DEVICE="${DEVICE:-cuda}"

# Ensure observed fasta present (copy from data if needed)
if [[ ! -f "$CASE_DIR/observed_anc_aa.fasta" ]]; then
  echo "ERROR: missing $CASE_DIR/observed_anc_aa.fasta — sync cases or run TreeSBM job first"
  exit 1
fi

$PYTHON scripts/gen_voc_baseline_leaves.py \
  --case-dir "$CASE_DIR" \
  --voc "$VOC_ID" \
  --methods $METHODS \
  --n-leaves "$N_LEAVES" \
  --device "$DEVICE"

echo "Done baselines for $VOC_ID @ $CASE_DIR"
