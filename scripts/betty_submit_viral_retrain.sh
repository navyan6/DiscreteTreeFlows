#!/bin/bash
# Submit viral mut+lineage retrains (new ckpt dirs; never overwrite paper ckpts).
# Recipe: v5/lit-hotspot family + β ∈ {0, 0.25} site_local.
#
#   covid_v6_mutlin[_b0|_b025]
#   h1n1_v3_mutlin[_b0|_b025]
#   h3n2_v4_mutlin[_b0|_b025]
#   hiv_geo_v2_mutlin[_b0|_b025]
#
# Usage: bash scripts/betty_submit_viral_retrain.sh

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs benchmarks/results/tables

IDS_JSON=benchmarks/results/tables/viral_retrain_mutlin_job_ids.json
echo "{" > "$IDS_JSON.tmp"
first=1

submit() {
  local name="$1"; shift
  local jid
  jid=$(sbatch --parsable "$@")
  echo "submitted $name -> $jid"
  if [[ $first -eq 1 ]]; then first=0; else echo "," >> "$IDS_JSON.tmp"; fi
  printf '  "%s": %s' "$name" "$jid" >> "$IDS_JSON.tmp"
}

# Guard: refuse if target would overwrite locked paper dirs
for d in covid_v5_mutrec h1n1_v2_lit_hotspot h3n2_v3_lit_hotspot hiv_geo_v1 ab_oas_1m_v1; do
  : # paper dirs are never used as CKPT_DIR below
done

for BETA in 0 0.25; do
  TAG="b0"
  [[ "$BETA" == "0.25" ]] && TAG="b025"

  # COVID Brazil geo + PMC lit mask (v5 recipe + optional β)
  submit "covid_v6_mutlin_${TAG}" \
    --qos=mig-max --job-name="cv6_${TAG}" --time=72:00:00 \
    --export=ALL,HOTSPOT=1,CKPT_DIR="checkpoints/covid_v6_mutlin_${TAG}",FITNESS_BETA="$BETA",FITNESS_SCORE=log_R0,MUT_HOTSPOT_MASK=results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt,MUT_HOTSPOT_FRAC=,MUT_HOTSPOT_TOPK=,MUT_HOTSPOT_WEIGHT=5,MUT_HOTSPOT_FORCE=1 \
    scripts/slurm_covid_train_mut_recovery.sh

  # H1N1 geo lit
  submit "h1n1_v3_mutlin_${TAG}" \
    --qos=mig-max --job-name="h1v3_${TAG}" --time=72:00:00 \
    --export=ALL,CKPT_DIR="checkpoints/h1n1_v3_mutlin_${TAG}",FITNESS_BETA="$BETA",FITNESS_SCORE=log_R0,MUT_HOTSPOT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_h1_nmicrobiol_lit.pt,MUT_HOTSPOT_WEIGHT=5,MUT_HOTSPOT_FORCE=1 \
    scripts/slurm_h1n1_train_lit_hotspot.sh

  # H3N2 temporal lit
  submit "h3n2_v4_mutlin_${TAG}" \
    --qos=mig-max --job-name="h3v4_${TAG}" --time=72:00:00 \
    --export=ALL,CKPT_DIR="checkpoints/h3n2_v4_mutlin_${TAG}",FITNESS_BETA="$BETA",FITNESS_SCORE=log_R0,MUT_HOTSPOT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt,MUT_HOTSPOT_WEIGHT=5,MUT_HOTSPOT_FORCE=1 \
    scripts/slurm_h3n2_train_lit_hotspot.sh

  # HIV geo V1–V5
  submit "hiv_geo_v2_mutlin_${TAG}" \
    --qos=mig-max --partition=b200-mig90 --job-name="hiv2_${TAG}" --time=72:00:00 \
    --export=ALL,DATA=data/hiv_geo,CKPT_DIR="checkpoints/hiv_geo_v2_mutlin_${TAG}",FITNESS_BETA="$BETA",FITNESS_SCORE=log_R0,SKIP_PRECOMPUTE=1 \
    scripts/slurm_hiv_train.sh
done

echo "" >> "$IDS_JSON.tmp"
echo "}" >> "$IDS_JSON.tmp"
mv "$IDS_JSON.tmp" "$IDS_JSON"
echo "Wrote $IDS_JSON"
cat "$IDS_JSON"
echo "NOTE: do not overwrite covid_v5_mutrec / h1n1_v2_lit_hotspot / h3n2_v3_lit_hotspot / hiv_geo_v1"
