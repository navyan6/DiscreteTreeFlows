#!/usr/bin/env bash
# Frame-audit one representative split per dataset family, so the translation
# mode for each is chosen from evidence rather than from the name.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
PY=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
D=/vast/projects/pranam/lab/nnori/DiscreteTreeFlows_data

DIRS=(
  data/bdbv_temporal/train
  data/bdbv_pan_temporal/train
  data/covid_voc_roots/uk_alpha
  "$D/covid_cladeholdout/train"
  "$D/covid_epidemic/train"
  "$D/covid_temporal/train"
  "$D/h1n1_epidemic/train"
  "$D/h1n1_temporal/train"
  "$D/h3n2_epidemic/train"
  "$D/h3n2_forecast_2020_cal/train"
  "$D/panflu_forecast_h3n2_cal/train"
  "$D/panflu_forecast_flub_cal/train"
  "$D/hiv_geo/val"
)

$PY -W ignore scripts/panviral/audit_frame_offset.py "${DIRS[@]}" --max-trees "${1:-10}" 2>&1 \
  | grep -E '===|frame call|protein median|gap runs|offset [0-9]:'
