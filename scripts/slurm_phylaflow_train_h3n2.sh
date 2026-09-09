#!/bin/bash
#SBATCH --job-name=pf_h3n2
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/vast/home/n/nnori/DiscreteTreeFlows/logs/phylaflow_train_h3n2_%j.log
#SBATCH --error=/vast/home/n/nnori/DiscreteTreeFlows/logs/phylaflow_train_h3n2_%j.log
#
# Train native PhylaFlow on H3N2 size-N bank (Table 2 — NOT DS1–8).
# Forbidden: ./launch_ds_local.sh ds1…ds8
#
# Usage:
#   sbatch --qos=mig-max scripts/slurm_phylaflow_train_h3n2.sh 16
#   sbatch --qos=mig-max scripts/slurm_phylaflow_train_h3n2.sh 32
#
# Prereqs:
#   1) Clone + venv at $LABHOME/baselines/PhylaFlow
#   2) Bank built: $PHYLAFLOW_DATA_ROOT/h3n2_N{N}/ + configs/h3n2_N{N}.yaml
#      (python benchmarks/external_adapters/build_h3n2_phylaflow_bank.py ...)

set -euo pipefail
N="${1:?usage: sbatch slurm_phylaflow_train_h3n2.sh <N>}"
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export DTF="${DTF:-/vast/home/n/nnori/DiscreteTreeFlows}"
export REPO_DIR="${LABHOME}/baselines/PhylaFlow"
export PHYLAFLOW_DATA_ROOT="${PHYLAFLOW_DATA_ROOT:-$LABHOME/baselines/phylaflow_h3n2_data}"
export PHYLAFLOW_ARTIFACT_ROOT="${PHYLAFLOW_ARTIFACT_ROOT:-$LABHOME/baselines/phylaflow_h3n2_artifacts}"
export PHYLAFLOW_OUTPUT_ROOT="${PHYLAFLOW_OUTPUT_ROOT:-$REPO_DIR/outputs_h3n2}"
mkdir -p "$DTF/logs" "$PHYLAFLOW_OUTPUT_ROOT" "$PHYLAFLOW_ARTIFACT_ROOT"

echo "=== PhylaFlow H3N2 train N=${N} (NOT DS1-8) ==="
echo "date=$(date) host=$(hostname)"
echo "REPO=$REPO_DIR"
echo "DATA_ROOT=$PHYLAFLOW_DATA_ROOT"
echo "OUTPUT_ROOT=$PHYLAFLOW_OUTPUT_ROOT"

if [ ! -d "$REPO_DIR/.git" ]; then
  echo "ERROR: PhylaFlow not cloned at $REPO_DIR"
  exit 1
fi
if [ ! -d "${PHYLAFLOW_DATA_ROOT}/h3n2_N${N}" ]; then
  echo "ERROR: missing H3N2 bank ${PHYLAFLOW_DATA_ROOT}/h3n2_N${N}"
  echo "Build with benchmarks/external_adapters/build_h3n2_phylaflow_bank.py"
  exit 1
fi
CFG="${REPO_DIR}/configs/h3n2_N${N}.yaml"
if [ ! -f "$CFG" ]; then
  echo "ERROR: missing config $CFG"
  exit 1
fi

# Activate PhylaFlow venv (created by setup job)
if [ -f "${REPO_DIR}/.venv/bin/activate" ]; then
  # shellcheck disable=SC1091
  source "${REPO_DIR}/.venv/bin/activate"
  PYTHON="${REPO_DIR}/.venv/bin/python"
elif [ -x /vast/home/n/nnori/.conda/envs/treesbm/bin/python ]; then
  PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
  echo "WARNING: using treesbm python — prefer PhylaFlow .venv"
else
  PYTHON=python
fi

export WANDB_MODE="${WANDB_MODE:-offline}"
cd "$REPO_DIR"
echo "python=$PYTHON"
$PYTHON -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"

echo "=== Launch train: python -m run.run $CFG ==="
# Explicit H3N2 config — never launch_ds_local.sh
$PYTHON -m run.run "$CFG"

echo "Done train N=${N}: $(date)"
echo "Next: sample dumps with evaluate_per_dataset_sample_kl.py --dump-trees"
echo "Then: sbatch --qos=mig-max $DTF/scripts/slurm_phylaflow.sh ${N}"
