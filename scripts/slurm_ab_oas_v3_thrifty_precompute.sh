#!/bin/bash
#SBATCH --job-name=ab_oas_v3_pre
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ab_oas_v3_thrifty_pre_%j.log
#SBATCH --error=logs/ab_oas_v3_thrifty_pre_%j.log
#
# Recipe B: Thrifty AA-marginal Q0 caches → group_*_ref_rates_thrifty.pt
# Does NOT overwrite ESM group_*_ref_rates.pt or PLM caches.
#
#   sbatch scripts/slurm_ab_oas_v3_thrifty_precompute.sh

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows
export PYTHONPATH="$PWD:${PWD}/antibody_benchmark/data/raw/repos/netam:${PYTHONPATH:-}"
mkdir -p logs

DATA="${DATA:-data/ab_clones_1m}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-160}"

if [[ ! -d "$DATA/train" ]]; then
  echo "ERROR: missing $DATA/train" >&2
  exit 1
fi

echo "=== Recipe B precompute thrifty_aa Q0 ==="
echo "data=$DATA max_seq_len=$MAX_SEQ_LEN host=$(hostname) date=$(date -Is)"

for split in train val test; do
  if [[ ! -d "$DATA/$split" ]]; then
    echo "skip missing $DATA/$split"
    continue
  fi
  echo "--- $split ---"
  $PYTHON -u scripts/precompute_ref_rates.py \
    --data "$DATA/$split" \
    --max-seq-len "$MAX_SEQ_LEN" \
    --r0-backend thrifty_aa
done

echo "Done: $(date -Is)"
ls -1 "$DATA/train"/group_*_ref_rates_thrifty.pt 2>/dev/null | wc -l
