#!/usr/bin/env bash
# Retranslate every dataset with the mode its own frame audit selected.
#
# Pass --dry-run first. covid/test is already repaired and must report zero
# changes; if it does not, the translator has drifted and nothing else should be
# written until that is understood.
#
#   bash scripts/panviral/retranslate_all.sh --dry-run
#   bash scripts/panviral/retranslate_all.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
PY=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
D=/vast/projects/pranam/lab/nnori/DiscreteTreeFlows_data

DIRS=(
  data/covid/train data/covid/val data/covid/test
  data/covid_voc_roots/*/
  data/bdbv_temporal/train data/bdbv_temporal/test
  data/bdbv_pan_temporal/train data/bdbv_pan_temporal/test
  "$D"/covid_cladeholdout/{train,val,test}
  "$D"/covid_epidemic/{train,val,test}
  "$D"/covid_temporal/{train,val,test}
  "$D"/h1n1/{train,val,test}
  "$D"/h1n1_epidemic/{train,val,test}
  "$D"/h1n1_temporal/{train,val,test}
  "$D"/h3n2_epidemic/{train,val,test}
  "$D"/h3n2_forecast_2020_cal/{train,test}
  "$D"/h3n2_forecast_2020_season/{train,test}
  "$D"/panflu_forecast_h3n2_cal/{train,test}
  "$D"/panflu_forecast_h3n2_cal_dual/{train,test}
  "$D"/panflu_forecast_flub_cal/train
  "$D"/hiv_geo/{train,val,test}
  "$D"/hiv_temporal/{train,val,test}
)

EXIST=()
for d in "${DIRS[@]}"; do
  [ -d "$d" ] && compgen -G "$d/group_*_anc_nt.fasta" > /dev/null && EXIST+=("$d")
done

echo "retranslating ${#EXIST[@]} dataset dirs ${1:-(WRITING)}"
$PY -W ignore scripts/retranslate_anc_aa.py --data-dir "${EXIST[@]}" "$@"
