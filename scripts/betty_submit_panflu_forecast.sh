#!/bin/bash
# Sync scripts + queue pan-flu forecast matrix on Betty (prep → pipeline → train).
set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
cd "${ROOT:-$HOME/DiscreteTreeFlows}"

echo "=== Pan-flu split prep (FluB ingest + holdout dirs) ==="
bash scripts/betty_prep_panflu_forecast.sh

submit_row() {
  local name="$1" data_root="$2" ckpt="$3" maxlen="${4:-566}"
  local j_pipe j_pre j_train
  j_pipe=$(sbatch --parsable --export=ALL,DATA_ROOT="$data_root" \
    scripts/slurm_epidemic_pipeline.sh)
  j_pre=$(sbatch --parsable --dependency=afterok:"$j_pipe" \
    --export=ALL,DATA_ROOT="$data_root",MAX_SEQ_LEN="$maxlen" \
    scripts/slurm_epidemic_precompute.sh)
  j_train=$(sbatch --parsable --dependency=afterok:"$j_pre" --qos=mig-max \
    --job-name="pf_${name}" \
    --export=ALL,DATA_ROOT="$data_root",CKPT_DIR="$ckpt",MAX_SEQ_LEN="$maxlen,HOTSPOT=1" \
    scripts/slurm_epidemic_train_mut_recovery.sh)
  echo "$name pipe=$j_pipe pre=$j_pre train=$j_train -> $ckpt"
}

# Start with B0 baseline; extend to full matrix once B0 pipeline validates
submit_row B0_h3n2_cal data/h3n2_forecast_2020_cal checkpoints/h3n2_only_forecast_2020_cal
submit_row B1_h3n2_season data/h3n2_forecast_2020_season checkpoints/h3n2_only_forecast_2020_season
submit_row G1_dual_h3n2_cal data/panflu_forecast_h3n2_cal_dual checkpoints/dual_forecast_h3n2_cal
submit_row G3_pan_h3n2_cal data/panflu_forecast_h3n2_cal checkpoints/panflu_forecast_h3n2_cal
submit_row G6_pan_flub_cal data/panflu_forecast_flub_cal checkpoints/panflu_forecast_flub_cal

echo "Pan-flu jobs submitted (subset first). Add G2/G4/G5 after B0 smoke OK."
