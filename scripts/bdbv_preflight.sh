#!/bin/bash
# Verify BDBV track artifacts exist before sbatch. Run on Betty (or locally after rsync).
set -euo pipefail

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
PYTHON="${PYTHON:-python3}"
[[ -x /vast/home/n/nnori/.conda/envs/treesbm/bin/python ]] && \
  PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

fail=0
need() {
  if [[ ! -e "$1" ]]; then
    echo "MISSING: $1" >&2
    fail=1
  else
    echo "OK: $1"
  fi
}

need data/bdbv/manifest.json
need data/bdbv/l_nt/all.fasta
need data/bdbv/l_window/all.fasta
need results/bdbv_l_conservation/window_config.json
need results/bdbv_l_mask/mut_hotspot_mask_lit.pt
need data/filo_l/train/filo_train_group_001.fasta
need data/filo_l/test/filo_test_group_001.fasta

if [[ -f data/bdbv/FILO_SPLIT_AUDIT.json ]]; then
  $PYTHON scripts/audit_filo_splits.py --data data/filo_l || true
fi

for s in scripts/slurm_bdbv_pipeline.sh scripts/slurm_bdbv_precompute.sh \
         scripts/slurm_bdbv_train_mut_recovery.sh scripts/slurm_bdbv_coverage.sh \
         scripts/slurm_eval_bdbv.sh scripts/betty_submit_bdbv.sh; do
  need "$s"
done

if [[ $fail -ne 0 ]]; then
  echo "Preflight FAILED — run: bash scripts/betty_resume_bdbv_ingest.sh" >&2
  exit 1
fi
echo "Preflight passed."
