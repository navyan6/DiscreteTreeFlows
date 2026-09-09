#!/bin/bash
# Enrichment eval on mutlin retrains (does not overwrite paper-ckpt evals).
# Usage on Betty: bash scripts/betty_submit_mutlin_enrich.sh

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs checkpoints/viral_eval_sweep benchmarks/results/tables

IDS=benchmarks/results/tables/viral_mutlin_enrich_job_ids.json
echo "{" > "$IDS.tmp"
first=1

submit() {
  local name="$1"; shift
  local ckpt="$1"; shift
  if [[ ! -f "$ckpt" ]]; then
    echo "SKIP $name (missing $ckpt)" >&2
    return 0
  fi
  local jid
  jid=$(sbatch --parsable "$@")
  echo "submitted $name -> $jid"
  if [[ $first -eq 1 ]]; then first=0; else echo "," >> "$IDS.tmp"; fi
  printf '  "%s": %s' "$name" "$jid" >> "$IDS.tmp"
}

MRS="${MRS:-0.5}"

# COVID
for TAG in b0 b025; do
  submit "covid_v6_mutlin_${TAG}_mrs${MRS}" \
    "checkpoints/covid_v6_mutlin_${TAG}/best.pt" \
    --qos=mig-max --job-name="en_cv6_${TAG}" \
    --export=ALL,CKPT="checkpoints/covid_v6_mutlin_${TAG}/best.pt",DATA=data/covid/test,MAX_SEQ_LEN=1280,MRS="$MRS",MAX_TREES=20,NSTEPS=100,OUT="checkpoints/viral_eval_sweep/covid_v6_mutlin_${TAG}_mrs${MRS}_enrich.json",LIT_MASK=results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt \
    scripts/slurm_c1_treesbm_plmnll.sh
done

# H1N1
for TAG in b0 b025; do
  submit "h1n1_v3_mutlin_${TAG}_mrs${MRS}" \
    "checkpoints/h1n1_v3_mutlin_${TAG}/best.pt" \
    --qos=mig-max --job-name="en_h1v3_${TAG}" \
    --export=ALL,CKPT="checkpoints/h1n1_v3_mutlin_${TAG}/best.pt",DATA=data/h1n1/test,MAX_SEQ_LEN=566,MRS="$MRS",MAX_TREES=20,NSTEPS=100,OUT="checkpoints/viral_eval_sweep/h1n1_v3_mutlin_${TAG}_mrs${MRS}_enrich.json",LIT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_h1_nmicrobiol_lit.pt \
    scripts/slurm_c1_treesbm_plmnll.sh
done

# H3N2
for TAG in b0 b025; do
  submit "h3n2_v4_mutlin_${TAG}_mrs${MRS}" \
    "checkpoints/h3n2_v4_mutlin_${TAG}/best.pt" \
    --qos=mig-max --job-name="en_h3v4_${TAG}" \
    --export=ALL,CKPT="checkpoints/h3n2_v4_mutlin_${TAG}/best.pt",DATA=data/h3n2/test,MAX_SEQ_LEN=566,MRS="$MRS",MAX_TREES=20,NSTEPS=100,OUT="checkpoints/viral_eval_sweep/h3n2_v4_mutlin_${TAG}_mrs${MRS}_enrich.json",LIT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt \
    scripts/slurm_c1_treesbm_plmnll.sh
done

# HIV geo
for TAG in b0 b025; do
  submit "hiv_geo_v2_mutlin_${TAG}_mrs${MRS}" \
    "checkpoints/hiv_geo_v2_mutlin_${TAG}/best.pt" \
    --qos=mig-max --job-name="en_hiv2_${TAG}" \
    --export=ALL,CKPT="checkpoints/hiv_geo_v2_mutlin_${TAG}/best.pt",DATA=data/hiv_geo/test,MAX_SEQ_LEN=900,MRS="$MRS",MAX_TREES=14,NSTEPS=100,OUT="checkpoints/viral_eval_sweep/hiv_geo_v2_mutlin_${TAG}_mrs${MRS}_enrich.json",LIT_MASK=results/hiv_env_mask/mut_hotspot_mask_hiv_env_v1v5.pt \
    scripts/slurm_c1_treesbm_plmnll.sh
done

echo "" >> "$IDS.tmp"
echo "}" >> "$IDS.tmp"
mv "$IDS.tmp" "$IDS"
echo "Wrote $IDS"
cat "$IDS"
squeue -u "$USER" | head -25
