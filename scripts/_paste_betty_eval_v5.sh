#!/bin/bash
# Paste on Betty (cwd ~/DiscreteTreeFlows) or: bash scripts/_paste_betty_eval_v5.sh
set -euo pipefail
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints/covid_v5_mutrec

CKPT=checkpoints/covid_v5_mutrec/best.pt
if [[ ! -f "$CKPT" && -f "$LABHOME/checkpoints_backup/covid_v5_mutrec/best.pt" ]]; then
  cp -n "$LABHOME/checkpoints_backup/covid_v5_mutrec/best.pt" "$CKPT"
fi
[[ -f "$CKPT" ]] || { echo "ERROR: missing $CKPT"; exit 1; }
ls -lah "$CKPT"

OUT03=checkpoints/eval_enrichment_covid_v5_mutrec_mrs0.3.json
OUT05=checkpoints/eval_enrichment_covid_v5_mutrec_mrs0.5.json
if [[ -f "$OUT03" && -f "$OUT05" ]]; then
  echo "SKIP: eval JSONs already exist"; echo "JOB_ID="; exit 0
fi
if squeue -u nnori -h -n eval_covid_v5 2>/dev/null | grep -q .; then
  echo "SKIP: eval_covid_v5 already queued/running:"
  squeue -u nnori -n eval_covid_v5
  echo "JOB_ID="; exit 0
fi

JOB=$(CKPT_NAME=covid_v5_mutrec EVAL_MRS="0.3 0.5" SKIP_EXISTING=1 \
  sbatch --qos=mig-max --parsable \
  --job-name=eval_covid_v5 \
  --output=logs/eval_covid_v5_%j.log \
  --error=logs/eval_covid_v5_%j.log \
  scripts/slurm_eval_covid_v4.sh)
echo "JOB_ID=$JOB"
