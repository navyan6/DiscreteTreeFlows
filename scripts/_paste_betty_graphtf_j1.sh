#!/bin/bash
# Paste on Betty: full GraphTF J.1 probe (NOT smoke). Writes non-smoke out dirs.
set -euo pipefail
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
cd ~/DiscreteTreeFlows
mkdir -p logs results

resolve_ckpt() {
  local name="$1"
  for p in \
    "checkpoints/${name}/best.pt" \
    "${LABHOME}/checkpoints_backup/${name}/best.pt"
  do
    if [[ -f "$p" ]]; then echo "$p"; return 0; fi
  done
  return 1
}

echo "=== GraphTF J.1 submit ==="
H3=$(resolve_ckpt h3n2_v2 || true)
COVID=$(resolve_ckpt covid_v5_mutrec || true)
echo "h3n2: ${H3:-MISSING}"
echo "covid_v5: ${COVID:-MISSING}"

QOS="${QOS:-mig-max}"

if [[ -n "${H3}" ]]; then
  if squeue -u nnori -h -n tf_val_h3n2 2>/dev/null | grep -q .; then
    echo "SKIP tf_val_h3n2 already queued"
  else
    # Dedicated dir — do not touch results/transformer_val_smoke
    J=$(sbatch --qos="$QOS" --parsable \
      --job-name=tf_val_h3n2 \
      --output=logs/tf_val_h3n2_%j.log \
      --error=logs/tf_val_h3n2_%j.log \
      scripts/slurm_validate_transformer.sh \
      "$H3" data/h3n2/train results/transformer_val_h3n2_j1 566)
    echo "tf_val_h3n2_JOB=$J"
  fi
else
  echo "ERROR: no h3n2_v2 ckpt" >&2
fi

# Optional COVID J.1 (longer seq); skip if data missing
if [[ -n "${COVID}" && -d data/covid/train ]]; then
  if squeue -u nnori -h -n tf_val_covid 2>/dev/null | grep -q .; then
    echo "SKIP tf_val_covid already queued"
  else
    J=$(sbatch --qos="$QOS" --parsable \
      --job-name=tf_val_covid \
      --output=logs/tf_val_covid_%j.log \
      --error=logs/tf_val_covid_%j.log \
      scripts/slurm_validate_transformer.sh \
      "$COVID" data/covid/train results/transformer_val_covid_v5_j1 1280)
    echo "tf_val_covid_v5_JOB=$J"
  fi
else
  echo "SKIP covid J.1 (ckpt or data/covid/train missing)"
fi

squeue -u nnori -o "%.18i %.12P %.10q %.20j %.8T %.10M %R" | head -40
echo "Primary tables: results/transformer_val_*_j1/probe_table_j1_primary.csv"
