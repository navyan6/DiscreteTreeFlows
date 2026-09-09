#!/bin/bash
# Prepare pan-flu forecast holdout matrix (B0–G7) on Betty.
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
cd "${ROOT:-$HOME/DiscreteTreeFlows}"

LAB_DATA="${LAB_DATA:-/vast/projects/pranam/lab/nnori/DiscreteTreeFlows_data}"
link_pool() {
  local name="$1"
  local repo_dir="data/${name}"
  local lab_dir="${LAB_DATA}/${name}"
  mkdir -p "$repo_dir"
  for split in train val test; do
    local target="${lab_dir}/${split}"
    local link="${repo_dir}/${split}"
    if [[ -d "$target" && ! -e "$link" ]]; then
      ln -s "$target" "$link"
      echo "linked ${link} -> ${target}"
    fi
  done
}
echo "=== Link H3N2/H1N1 pools from lab storage ==="
link_pool h3n2
link_pool h1n1

echo "=== Flu B NCBI ingest ==="
$PYTHON scripts/download_flub_ha_ncbi.py --retmax 6000 --max-records 8000

echo "=== Pan-flu pool ==="
$PYTHON scripts/prepare_panflu_pool.py --include-flub

prep() {
  python scripts/prepare_panflu_forecast.py "$@"
}

# Row B0 — H3N2-only → calendar 2020 (baseline)
prep --train-pool h3n2 --test-subtype h3n2 --test-mode calendar \
  --out-base data/h3n2_forecast_2020_cal --alias-out data/h3n2_forecast_2020_cal

# Row B1 — H3N2-only → NH season 2019-2020
prep --train-pool h3n2 --test-subtype h3n2 --test-mode season \
  --out-base data/h3n2_forecast_2020_season

# G1/G2 — dual strain
prep --train-pool h3n2_h1n1 --test-subtype h3n2 --test-mode calendar \
  --out-base data/panflu_forecast_h3n2_cal_dual
prep --train-pool h3n2_h1n1 --test-subtype h1n1 --test-mode calendar \
  --out-base data/panflu_forecast_h1n1_cal_dual

# G3–G6 — pan-flu train
prep --train-pool panflu --test-subtype h3n2 --test-mode calendar \
  --out-base data/panflu_forecast_h3n2_cal
prep --train-pool panflu --test-subtype h3n2 --test-mode season \
  --out-base data/panflu_forecast_h3n2_season
prep --train-pool panflu --test-subtype h1n1 --test-mode calendar \
  --out-base data/panflu_forecast_h1n1_cal
prep --train-pool panflu --test-subtype flub --test-mode calendar \
  --out-base data/panflu_forecast_flub_cal

python scripts/audit_epidemic_splits.py --dirs \
  data/h3n2_forecast_2020_cal data/h3n2_forecast_2020_season \
  data/panflu_forecast_h3n2_cal data/panflu_forecast_h3n2_season || true

echo "Pan-flu forecast dirs ready under data/panflu_forecast_* and data/h3n2_forecast_*"
