#!/bin/bash
# One-shot: sync already done locally; run ON cluster after scp, or run via:
#   bash scripts/_cluster_submit_graphtf_probe.sh
set -euo pipefail
cd "${HOME}/DiscreteTreeFlows"
LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
mkdir -p logs results
chmod +x scripts/validate_transformer.py scripts/slurm_validate_transformer.sh \
         scripts/check_column_entropy.py scripts/slurm_check_column_entropy.sh

resolve_ckpt() {
  local name="$1"
  for p in \
    "checkpoints/${name}/best.pt" \
    "checkpoints/${name}/best_model.pt" \
    "${LABHOME}/checkpoints_backup/${name}/best.pt" \
    "${LABHOME}/checkpoints_backup/${name}/best_model.pt" \
    "${LABHOME}/checkpoints_backup/${name}/checkpoint_best.pt"
  do
    if [[ -f "$p" ]]; then echo "$p"; return 0; fi
  done
  return 1
}

echo "=== submitting GraphTF validate + optional entropy ==="
H3=$(resolve_ckpt h3n2_v2 || true)
COVID=$(resolve_ckpt covid_v3_cons || true)
echo "h3n2 ckpt: ${H3:-MISSING}"
echo "covid ckpt: ${COVID:-MISSING}"

QOS="${QOS:-mig-max}"
echo "Using QOS=${QOS}"

if [[ -n "${H3}" ]]; then
  sbatch --qos="${QOS}" scripts/slurm_validate_transformer.sh \
    "$H3" data/h3n2/train results/transformer_val_h3n2_j1 566
else
  echo "ERROR: no h3n2_v2 checkpoint found" >&2
fi

if [[ -n "${COVID}" && -d data/covid/train ]]; then
  sbatch --qos="${QOS}" scripts/slurm_validate_transformer.sh \
    "$COVID" data/covid/train results/transformer_val_covid 1280
else
  echo "SKIP covid validate (ckpt or data/covid/train missing)"
fi

if [[ -d data/covid/train ]]; then
  sbatch scripts/slurm_check_column_entropy.sh \
    data/covid/train 1280 results/column_entropy_covid
else
  echo "SKIP covid entropy (data/covid/train missing)"
fi

echo "=== squeue ==="
squeue -u "$USER" -o '%.18i %.12P %.10q %.16j %.8T %.10M %R' | head -40
