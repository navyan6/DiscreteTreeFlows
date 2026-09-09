#!/bin/bash
# Full antibody benchmark rollout (one model per job).
#
# CPU models (thrifty / dasm_thrifty / treesbm):
#   MODEL=thrifty sbatch antibody_benchmark/scripts/slurm_ab_full_rollout.sh
#
# CoSiNE GPU (prefer mig-max on B200 MIG):
#   MODEL=cosine sbatch --partition=b200-mig45 --qos=mig-max --gres=gpu:1 \
#     --cpus-per-task=6 --mem-per-cpu=8G \
#     antibody_benchmark/scripts/slurm_ab_full_rollout.sh
#
#SBATCH --job-name=ab_full
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=/vast/projects/pranam/lab/nnori/antibody_benchmark/logs/ab_full_%x_%j.out
#SBATCH --error=/vast/projects/pranam/lab/nnori/antibody_benchmark/logs/ab_full_%x_%j.err

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"
export TORCH_HOME="${TORCH_HOME:-$LABHOME/torch_cache}"

ROOT="${ROOT:-$HOME/DiscreteTreeFlows}"
cd "$ROOT"
MODEL="${MODEL:?Set MODEL=neutral_shm|plm_prior|ar_tree_edit|thrifty|dasm_thrifty|cosine|treesbm}"
CONFIG="${CONFIG:-antibody_benchmark/configs/full.yaml}"
PY="${PY:-/vast/projects/pranam/lab/nnori/.conda/envs/treesbm/bin/python}"
N_ROLLOUTS="${N_ROLLOUTS:-20}"

mkdir -p "$LABHOME/antibody_benchmark/logs" \
  antibody_benchmark/results/samples \
  antibody_benchmark/results/summary

export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
export AB_DASM_DIR="${AB_DASM_DIR:-$HOME/antibody_benchmark_raw/dasm/extracted/dasm-experiments-data}"
if [[ -d antibody_benchmark/data/raw/repos/cosine ]]; then
  export PYTHONPATH="$ROOT/antibody_benchmark/data/raw/repos/cosine:$ROOT/antibody_benchmark/data/raw/repos/cosine/evo:${PYTHONPATH}"
fi

# CoSiNE ckpt: prefer LABHOME canonical, symlink into repo path if needed
CKPT_REPO=antibody_benchmark/data/raw/cosine/checkpoints/cosine_dasm.ckpt
CKPT_LAB="$LABHOME/antibody_benchmark/cosine_ckpts/cosine_dasm.ckpt"
if [[ "$MODEL" == "cosine" || "$MODEL" == "plm_prior" ]]; then
  module load cuda/12.8 2>/dev/null || module load cuda 2>/dev/null || true
fi

if [[ "$MODEL" == "cosine" ]]; then
  mkdir -p antibody_benchmark/data/raw/cosine/checkpoints
  if [[ ! -f "$CKPT_REPO" && -f "$CKPT_LAB" ]]; then
    ln -sfn "$CKPT_LAB" "$CKPT_REPO"
  fi
  if [[ ! -f "$CKPT_REPO" ]]; then
    echo "BLOCKER: missing CoSiNE ckpt ($CKPT_REPO or $CKPT_LAB)"
    exit 1
  fi
fi

if [[ "$MODEL" == "ar_tree_edit" ]]; then
  POOL=benchmarks/external_pools/sampled/artreeformer_N16.nwk
  if [[ ! -f "$POOL" ]]; then
    echo "BLOCKER: missing ARTreeFormer pool $POOL"
    exit 1
  fi
fi

if [[ "$MODEL" == "treesbm" && ! -f checkpoints/best.pt ]]; then
  echo "BLOCKER: missing TreeSBM checkpoint checkpoints/best.pt"
  exit 1
fi

echo "=== ab full rollout ==="
echo "host=$(hostname) date=$(date -Is) model=$MODEL n_rollouts=$N_ROLLOUTS"
echo "root=$ROOT config=$CONFIG"
echo "cuda=${CUDA_VISIBLE_DEVICES:-none}"
"$PY" - <<'PY'
import torch
print("torch", torch.__version__, "cuda", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device", torch.cuda.get_device_name(0))
PY

"$PY" antibody_benchmark/scripts/run_rollouts.py \
  --config "$CONFIG" \
  --models "$MODEL" \
  --n-rollouts "$N_ROLLOUTS"

# Count samples written
nsamp=$(find "antibody_benchmark/results/samples/$MODEL" -name 'rollout_*.json' 2>/dev/null | wc -l | tr -d ' ')
echo "samples_written=$nsamp model=$MODEL"
echo "Done: $(date -Is)"
