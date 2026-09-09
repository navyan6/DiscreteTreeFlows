#!/bin/bash
#SBATCH --job-name=epidemic_pre
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/epidemic_precompute_%j.log
#SBATCH --error=logs/epidemic_precompute_%j.log
#
# Env: DATA_ROOT, MAX_SEQ_LEN (566 flu, 1280 covid, 900 filo)

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows
mkdir -p logs

DATA_ROOT="${DATA_ROOT:?set DATA_ROOT}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"
if [[ "$DATA_ROOT" == *filo* ]] && [[ -f results/bdbv_l_conservation/window_config.json ]]; then
  MAX_SEQ_LEN=$($PYTHON -c "import json; print(json.load(open('results/bdbv_l_conservation/window_config.json'))['max_seq_len'])")
fi

echo "precompute DATA_ROOT=$DATA_ROOT MAX_SEQ_LEN=$MAX_SEQ_LEN"

for split in train val test; do
  d="$DATA_ROOT/$split"
  [[ -d "$d" ]] || continue
  [[ -n "$(ls -A "$d"/*_group_*.fasta 2>/dev/null || true)" ]] || continue
  echo "=== PLM $split ==="
  $PYTHON scripts/precompute_plm.py --data "$d"
  echo "=== ref_rates $split ==="
  $PYTHON scripts/precompute_ref_rates.py --data "$d" --max-seq-len "$MAX_SEQ_LEN"
done

echo "Done: $(date)"
