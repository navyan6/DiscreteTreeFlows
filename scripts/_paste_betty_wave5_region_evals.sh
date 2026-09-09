#!/bin/bash
# Paste on Betty: full-metric enrichment evals (region site_recall/any_mut/mut_frac)
# for latest Wave-5 trains at mrs ∈ {0.3, 0.5}. Writes *_regions.json so old
# JSONs (without region_site_recall) are not SKIP_EXISTING'd incorrectly.
set -euo pipefail
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

echo "=== BEFORE ==="
squeue -u nnori -o "%.18i %.12P %.20j %.8T %.10M %R %q" | head -40

submit_covid() {
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
  job=$(CKPT_NAME="$name" EVAL_MRS="0.3 0.5" OUT_SUFFIX=_regions SKIP_EXISTING=1 \
    sbatch --qos=mig-max --parsable \
      --job-name="$jname" \
      --output="logs/${jname}_%j.log" \
      --error="logs/${jname}_%j.log" \
      scripts/slurm_eval_covid_v4.sh)
  echo "${name}_JOB=$job"
}

submit_h3n2() {
  local name="h3n2_v3_lit_hotspot" jname="eval_h3n2_regions"
  local ckpt="checkpoints/${name}/best.pt"
  [[ -f "$ckpt" ]] || { echo "ERROR: missing $ckpt"; return 1; }
  ls -lah "$ckpt"
  if squeue -u nnori -h -n "$jname" 2>/dev/null | grep -q .; then
    echo "SKIP: $jname already queued"
    squeue -u nnori -n "$jname"
    echo "${name}_JOB="; return 0
  fi
  local job
  job=$(CKPT_NAME="$name" EVAL_MRS="0.3 0.5" OUT_SUFFIX=_regions SKIP_EXISTING=1 \
    sbatch --qos=mig-max --parsable \
      --job-name="$jname" \
      --output="logs/${jname}_%j.log" \
      --error="logs/${jname}_%j.log" \
      scripts/_run_eval_h3n2_lit.sh)
  echo "${name}_JOB=$job"
}

echo "=== SUBMIT Wave-5 region enrichments (OUT_SUFFIX=_regions) ==="
submit_covid covid_v7_pmc_hotspot eval_v7_regions
submit_covid covid_v5_deepmut eval_v5_deepmut_reg
submit_covid covid_v5_mutrec eval_v5_mutrec_reg
submit_h3n2

echo "=== AFTER ==="
squeue -u nnori -o "%.18i %.12P %.20j %.8T %.10M %R %q" | head -40
echo "Artifacts: checkpoints/eval_enrichment_<ckpt>_mrs{0.3,0.5}_regions.json"
echo "Note: covid_v7 mrs0.3_regions already exists (job 7437505) → SKIP_EXISTING will skip it."
