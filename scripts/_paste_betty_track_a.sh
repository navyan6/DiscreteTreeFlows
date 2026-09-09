#!/bin/bash
# Paste on Betty (cwd ~/DiscreteTreeFlows) or: bash scripts/_paste_betty_track_a.sh
# Submits Coverage@K track_a for covid_v5_mutrec (+ v6 if ckpt present).
set -euo pipefail
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results checkpoints/covid_v5_mutrec checkpoints/covid_v6_hotspot

resolve_ckpt() {
  local name="$1"
  if [[ -f "checkpoints/${name}/best.pt" ]]; then
    echo "checkpoints/${name}/best.pt"; return 0
  fi
  if [[ -f "$LABHOME/checkpoints_backup/${name}/best.pt" ]]; then
    mkdir -p "checkpoints/${name}"
    cp -n "$LABHOME/checkpoints_backup/${name}/best.pt" "checkpoints/${name}/best.pt"
    echo "checkpoints/${name}/best.pt"; return 0
  fi
  return 1
}

submit_one() {
  local name="$1" jname="$2"
  local ckpt
  ckpt=$(resolve_ckpt "$name") || { echo "SKIP $name (no ckpt)"; echo "${name}_JOB="; return 0; }
  ls -lah "$ckpt"
  if squeue -u nnori -h -n "$jname" 2>/dev/null | grep -q .; then
    echo "SKIP: $jname already queued/running:"
    squeue -u nnori -n "$jname"
    echo "${name}_JOB="; return 0
  fi
  local job
  job=$(CKPT_NAME="$name" DATA=data/covid/test MAX_SEQ_LEN=1280 KS="100 500 1000" \
    N_GROUPS=20 SKIP_EXISTING=1 \
    sbatch --qos=mig-max --parsable \
      --job-name="$jname" \
      --output="logs/${jname}_%j.log" \
      --error="logs/${jname}_%j.log" \
      scripts/slurm_track_a.sh)
  echo "${name}_JOB=$job"
}

echo "=== track_a submit $(date) ==="
submit_one covid_v5_mutrec track_a_v5
submit_one covid_v6_hotspot track_a_v6
echo "=== squeue ==="
squeue -u nnori -o "%.18i %.12P %.10q %.20j %.8T %.10M %R" | head -40
