#!/bin/bash
# Paste on Betty: force re-enrichment so JSONs include site_recall / aa_acc_given_hit.
# Cancels broken eval_v5_siteaa / eval_v6_siteaa if they lack CKPT_NAME / SKIP_EXISTING=0.
set -euo pipefail
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

echo "=== BEFORE ==="
squeue -u nnori -o "%.18i %.12P %.20j %.8T %.10M %R %q" | head -40

# Broken prior submits (no env on SubmitLine → default covid_v4 + SKIP_EXISTING=1)
for j in 7429564 7429565; do
  if squeue -j "$j" -h >/dev/null 2>&1; then
    echo "scancel $j (broken siteaa — missing CKPT_NAME/SKIP_EXISTING=0)"
    scancel "$j" || true
  fi
done
sleep 2

submit_reeval() {
  local name="$1" jname="$2"
  local ckpt="checkpoints/${name}/best.pt"
  if [[ ! -f "$ckpt" && -f "$LABHOME/checkpoints_backup/${name}/best.pt" ]]; then
    mkdir -p "checkpoints/${name}"
    cp -n "$LABHOME/checkpoints_backup/${name}/best.pt" "$ckpt"
  fi
  [[ -f "$ckpt" ]] || { echo "ERROR: missing $ckpt"; return 1; }
  ls -lah "$ckpt"
  if squeue -u nnori -h -n "$jname" 2>/dev/null | grep -q .; then
    echo "SKIP: $jname already queued"
    squeue -u nnori -n "$jname"
    echo "${name}_JOB="; return 0
  fi
  local job
  job=$(CKPT_NAME="$name" EVAL_MRS="0.3 0.5" SKIP_EXISTING=0 \
    sbatch --qos=mig-max --parsable \
      --job-name="$jname" \
      --output="logs/${jname}_%j.log" \
      --error="logs/${jname}_%j.log" \
      scripts/slurm_eval_covid_v4.sh)
  echo "${name}_JOB=$job"
}

echo "=== SUBMIT siteaa re-enrichment (overwrite JSONs) ==="
submit_reeval covid_v5_mutrec eval_v5_siteaa
submit_reeval covid_v6_hotspot eval_v6_siteaa

echo "=== AFTER ==="
squeue -u nnori -o "%.18i %.12P %.20j %.8T %.10M %R %q" | head -40
echo "Artifacts (after finish): checkpoints/eval_enrichment_<ckpt>_mrs{0.3,0.5}.json"
echo "Expect summary keys: site_recall, aa_acc_given_hit"
