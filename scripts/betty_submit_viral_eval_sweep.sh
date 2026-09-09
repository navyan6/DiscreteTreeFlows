#!/bin/bash
# Submit eval-only mrs / K / ε sweeps on current paper ckpts (no train).
# Primary metrics: mut_recovery, aa|hit, clade_recall, mean_min_edit.
# Cov@ε reported as secondary (ε∈{1,2,3}, K∈{100,500}).
#
# Usage (on Betty login03):
#   bash scripts/betty_submit_viral_eval_sweep.sh
#
# Writes job IDs to benchmarks/results/tables/viral_eval_sweep_job_ids.json

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs benchmarks/results/tables checkpoints/viral_eval_sweep

IDS_JSON=benchmarks/results/tables/viral_eval_sweep_job_ids.json
echo "{" > "$IDS_JSON.tmp"
first=1

submit() {
  local name="$1"; shift
  local jid
  jid=$(sbatch --parsable "$@")
  echo "submitted $name -> $jid"
  if [[ $first -eq 1 ]]; then
    first=0
  else
    echo "," >> "$IDS_JSON.tmp"
  fi
  printf '  "%s": %s' "$name" "$jid" >> "$IDS_JSON.tmp"
}

# --- Enrichment mrs grid (mut / aa|hit / clade) ---
# COVID
for MRS in 0.3 0.5 1.0; do
  submit "covid_enrich_mrs${MRS}" \
    --qos=mig-max --job-name="ve_cov_m${MRS}" \
    --export=ALL,CKPT=checkpoints/covid_v5_mutrec/best.pt,DATA=data/covid/test,MAX_SEQ_LEN=1280,MRS="$MRS",MAX_TREES=20,NSTEPS=100,OUT="checkpoints/viral_eval_sweep/covid_v5_mrs${MRS}_enrich.json",LIT_MASK=results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt \
    scripts/slurm_c1_treesbm_plmnll.sh
done

# H1N1
for MRS in 0.3 0.5 1.0; do
  submit "h1n1_enrich_mrs${MRS}" \
    --qos=mig-max --job-name="ve_h1_m${MRS}" \
    --export=ALL,CKPT=checkpoints/h1n1_v2_lit_hotspot/best.pt,DATA=data/h1n1/test,MAX_SEQ_LEN=566,MRS="$MRS",MAX_TREES=20,NSTEPS=100,OUT="checkpoints/viral_eval_sweep/h1n1_v2_mrs${MRS}_enrich.json",LIT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_h1_nmicrobiol_lit.pt \
    scripts/slurm_c1_treesbm_plmnll.sh
done

# H3N2
for MRS in 0.3 0.5 1.0; do
  submit "h3n2_enrich_mrs${MRS}" \
    --qos=mig-max --job-name="ve_h3_m${MRS}" \
    --export=ALL,CKPT=checkpoints/h3n2_v3_lit_hotspot/best.pt,DATA=data/h3n2/test,MAX_SEQ_LEN=566,MRS="$MRS",MAX_TREES=20,NSTEPS=100,OUT="checkpoints/viral_eval_sweep/h3n2_v3_mrs${MRS}_enrich.json",LIT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt \
    scripts/slurm_c1_treesbm_plmnll.sh
done

# HIV geo
for MRS in 0.3 0.5 1.0; do
  submit "hiv_geo_enrich_mrs${MRS}" \
    --qos=mig-max --job-name="ve_hiv_m${MRS}" \
    --export=ALL,CKPT=checkpoints/hiv_geo_v1/best.pt,DATA=data/hiv_geo/test,MAX_SEQ_LEN=900,MRS="$MRS",MAX_TREES=14,NSTEPS=100,OUT="checkpoints/viral_eval_sweep/hiv_geo_v1_mrs${MRS}_enrich.json",LIT_MASK=results/hiv_env_mask/mut_hotspot_mask_hiv_env_v1v5.pt \
    scripts/slurm_c1_treesbm_plmnll.sh
done

# --- Coverage K=100 and K=500 (ε abs Hamming; TreeSBM only) ---
for VIRUS in covid h1n1 hiv_geo; do
  case "$VIRUS" in
    covid) CK=checkpoints/covid_v5_mutrec/best.pt ;;
    h1n1) CK=checkpoints/h1n1_v2_lit_hotspot/best.pt ;;
    hiv_geo) CK=checkpoints/hiv_geo_v1/best.pt ;;
  esac
  for K in 100 500; do
    submit "${VIRUS}_cov_K${K}" \
      --qos=mig-max --job-name="ve_${VIRUS}_k${K}" \
      --export=ALL,VIRUS="$VIRUS",CHECKPOINT="$CK",METHODS=treesbm,K_MAX="$K",K_STEP=50,MAX_ROOTS=5,N_LIST=16,OUT="benchmarks/results/viral_eval_sweep_${VIRUS}_K${K}.csv",CACHE_DIR="benchmarks/results/coverage_cache_${VIRUS}_viral_eval_K${K}" \
      scripts/slurm_table5_coverage.sh
  done
done

# H3N2 coverage via coverage_curves script if table5 lacks h3n2 — use table5 with VIRUS override if supported
if grep -q 'h3n2' scripts/slurm_table5_coverage.sh; then
  for K in 100 500; do
    submit "h3n2_cov_K${K}" \
      --qos=mig-max --job-name="ve_h3_k${K}" \
      --export=ALL,VIRUS=h3n2,CHECKPOINT=checkpoints/h3n2_v3_lit_hotspot/best.pt,METHODS=treesbm,K_MAX="$K",K_STEP=50,MAX_ROOTS=5,N_LIST=16,OUT="benchmarks/results/viral_eval_sweep_h3n2_K${K}.csv",CACHE_DIR="benchmarks/results/coverage_cache_h3n2_viral_eval_K${K}" \
      scripts/slurm_table5_coverage.sh
  done
else
  # Fallback: reuse coverage_curves with CHECKPOINT env if script supports it
  submit "h3n2_cov_K100" \
    --qos=mig-max --job-name="ve_h3_k100" \
    --export=ALL,VIRUS=h3n2,CHECKPOINT=checkpoints/h3n2_v3_lit_hotspot/best.pt,METHODS=treesbm,K_MAX=100,OUT=benchmarks/results/viral_eval_sweep_h3n2_K100.csv \
    scripts/slurm_coverage_curves.sh || true
fi

echo "" >> "$IDS_JSON.tmp"
echo "}" >> "$IDS_JSON.tmp"
mv "$IDS_JSON.tmp" "$IDS_JSON"
echo "Wrote $IDS_JSON"
cat "$IDS_JSON"
