#!/bin/bash
#SBATCH --job-name=bdbv_pre
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/bdbv_precompute_%j.log
#SBATCH --error=logs/bdbv_precompute_%j.log
#
# MAX_SEQ_LEN from results/bdbv_l_conservation/window_config.json (default 900)

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows
mkdir -p logs

DATA_ROOT="${DATA_ROOT:-data/bdbv_temporal}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-900}"
if [[ -f results/bdbv_l_conservation/window_config.json ]]; then
  MAX_SEQ_LEN=$($PYTHON -c "import json; print(json.load(open('results/bdbv_l_conservation/window_config.json'))['max_seq_len'])")
fi

echo "precompute DATA_ROOT=$DATA_ROOT MAX_SEQ_LEN=$MAX_SEQ_LEN"

for split in train val test; do
  d="$DATA_ROOT/$split"
  [[ -d "$d" ]] || continue
  echo "=== PLM $split ==="
  $PYTHON scripts/precompute_plm.py --data "$d"
  echo "=== ref_rates $split ==="
  $PYTHON scripts/precompute_ref_rates.py --data "$d" --max-seq-len "$MAX_SEQ_LEN"
done

echo "Done: $(date)"
