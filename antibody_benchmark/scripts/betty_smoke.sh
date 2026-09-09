#!/usr/bin/env bash
# Run on Betty login/compute node after syncing antibody_benchmark/.
set -euo pipefail
ROOT="${ROOT:-$HOME/DiscreteTreeFlows}"
cd "$ROOT"
PY="${PY:-/vast/projects/pranam/lab/nnori/.conda/envs/treesbm/bin/python}"
export PYTHONPATH="$ROOT:${PYTHONPATH:-}"
export AB_DASM_DIR="${AB_DASM_DIR:-$HOME/antibody_benchmark_raw/dasm/extracted/dasm-experiments-data}"
# Optional: cosine repo + evo on path
if [[ -d antibody_benchmark/data/raw/repos/cosine ]]; then
  export PYTHONPATH="$ROOT/antibody_benchmark/data/raw/repos/cosine:$ROOT/antibody_benchmark/data/raw/repos/cosine/evo:${PYTHONPATH}"
fi
mkdir -p antibody_benchmark/results/smoke_reports
for m in thrifty dasm_thrifty cosine treesbm; do
  echo "=== SMOKE $m ==="
  "$PY" antibody_benchmark/scripts/smoke_models.py \
    --config antibody_benchmark/configs/default.yaml \
    --models "$m" --family-index 0 \
    || echo "FAILED $m"
done
echo "=== SUMMARY ==="
cat antibody_benchmark/results/smoke_reports/summary.json
