#!/bin/bash
#SBATCH --job-name=pf_setup
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=/vast/home/n/nnori/DiscreteTreeFlows/logs/phylaflow_setup_h3n2_%j.log
#SBATCH --error=/vast/home/n/nnori/DiscreteTreeFlows/logs/phylaflow_setup_h3n2_%j.log
#
# One-shot: PhylaFlow venv + H3N2 bank (N=16,32) for Table 2.
# Does NOT touch DS1–8 / launch_ds_local.sh.

set -euo pipefail
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export DTF="${DTF:-/vast/home/n/nnori/DiscreteTreeFlows}"
export REPO_DIR="${LABHOME}/baselines/PhylaFlow"
export PHYLAFLOW_DATA_ROOT="${PHYLAFLOW_DATA_ROOT:-$LABHOME/baselines/phylaflow_h3n2_data}"
export PHYLAFLOW_ARTIFACT_ROOT="${PHYLAFLOW_ARTIFACT_ROOT:-$LABHOME/baselines/phylaflow_h3n2_artifacts}"
export PHYLAFLOW_OUTPUT_ROOT="${PHYLAFLOW_OUTPUT_ROOT:-$REPO_DIR/outputs_h3n2}"
mkdir -p "$DTF/logs" "$PHYLAFLOW_DATA_ROOT" "$PHYLAFLOW_ARTIFACT_ROOT" "$PHYLAFLOW_OUTPUT_ROOT"

echo "=== PhylaFlow H3N2 setup ==="
echo "date=$(date) host=$(hostname)"
echo "REPO=$REPO_DIR DATA_ROOT=$PHYLAFLOW_DATA_ROOT"

if [ ! -d "$REPO_DIR/.git" ]; then
  git clone https://github.com/yashaektefaie/PhylaFlow "$REPO_DIR"
fi
(cd "$REPO_DIR" && git rev-parse HEAD | tee "$DTF/logs/phylaflow_clone_pin.txt")

cp -f "$DTF/benchmarks/external_adapters/phylaflow_sample.py" "$REPO_DIR/phylaflow_sample.py"
cp -f "$DTF/benchmarks/external_adapters/build_h3n2_phylaflow_bank.py" \
  "$REPO_DIR/build_h3n2_phylaflow_bank.py"

# Login image lacks ensurepip (python3 -m venv fails). Bootstrap with treesbm
# python + virtualenv, keep the env on LABHOME (home is 50GB).
TREESBM_PY="${TREESBM_PY:-/vast/home/n/nnori/.conda/envs/treesbm/bin/python}"
PHYLA_ENV="${PHYLAFLOW_VENV:-$LABHOME/baselines/conda_envs/phylaflow}"
mkdir -p "$(dirname "$PHYLA_ENV")"
if [ ! -x "${PHYLA_ENV}/bin/python" ]; then
  echo "=== Creating virtualenv at $PHYLA_ENV via treesbm python ==="
  "$TREESBM_PY" -m pip install --upgrade virtualenv
  "$TREESBM_PY" -m virtualenv "$PHYLA_ENV"
fi
# shellcheck disable=SC1091
source "${PHYLA_ENV}/bin/activate"
python -m pip install -U pip wheel
# Skip mamba-ssm (only for regenerating Phyla embeddings from raw seqs)
grep -v '^mamba-ssm' "${REPO_DIR}/requirements.txt" > /tmp/phylaflow_req_nomamba.txt
python -m pip install -r /tmp/phylaflow_req_nomamba.txt
python -m pip install "pytorch-lightning" "pyyaml" "ete3" || true
python -c "import torch, ete3, yaml; print('ok', torch.__version__)"
# Train script looks for $REPO_DIR/.venv
rm -rf "${REPO_DIR}/.venv"
ln -sfn "$PHYLA_ENV" "${REPO_DIR}/.venv"
echo "PHYLA_ENV=$PHYLA_ENV" | tee "$DTF/logs/phylaflow_env_path.txt"

POOL="$DTF/benchmarks/external_pools"
for N in 16 32; do
  TOPO="$POOL/train_topologies_N${N}.nwk"
  TRP="$POOL/train_topologies_N${N}.trprobs"
  if [ ! -f "$TOPO" ]; then
    echo "ERROR: missing $TOPO — run export_train_topologies.py first"
    exit 1
  fi
  echo "=== Building H3N2 bank N=${N} ==="
  python "$REPO_DIR/build_h3n2_phylaflow_bank.py" \
    --N "$N" \
    --topologies "$TOPO" \
    --trprobs "$TRP" \
    --n-cases 42 \
    --data-root "$PHYLAFLOW_DATA_ROOT" \
    --repo-dir "$REPO_DIR"
done

echo "=== Setup done $(date) ==="
ls -la "$PHYLAFLOW_DATA_ROOT"/h3n2_N16 "$PHYLAFLOW_DATA_ROOT"/h3n2_N32
ls -la "$REPO_DIR/configs"/h3n2_N*.yaml
echo "Next: sbatch --qos=mig-max $DTF/scripts/slurm_phylaflow_train_h3n2.sh 16"
echo "      sbatch --qos=mig-max $DTF/scripts/slurm_phylaflow_train_h3n2.sh 32"
