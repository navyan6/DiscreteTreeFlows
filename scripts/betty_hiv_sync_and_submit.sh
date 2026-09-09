#!/bin/bash
# Sync HIV prep artifacts to Betty and submit pipeline + trains.
# Run on a machine that can reach login.betty.parcc.upenn.edu (VPN).
#
# Usage (from repo root on laptop):
#   bash scripts/betty_hiv_sync_and_submit.sh
#
# Env:
#   BETTY=nnori@login.betty.parcc.upenn.edu
#   REMOTE=~/DiscreteTreeFlows
#   SKIP_RAW=1   # skip rsync of raw data/hiv/*.fasta if already on Betty
#   SUBMIT=1     # default 1; set 0 to sync only

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

BETTY="${BETTY:-nnori@login.betty.parcc.upenn.edu}"
REMOTE="${REMOTE:-~/DiscreteTreeFlows}"
SKIP_RAW="${SKIP_RAW:-0}"
SUBMIT="${SUBMIT:-1}"

echo "=== sync scripts ==="
rsync -avP \
  scripts/hiv_env_utils.py \
  scripts/prepare_hiv_splits.py \
  scripts/build_hiv_env_v1v5_hotspot_mask.py \
  scripts/slurm_hiv_pipeline.sh \
  scripts/slurm_hiv_train.sh \
  scripts/betty_hiv_sync_and_submit.sh \
  "$BETTY:$REMOTE/scripts/"

echo "=== sync mask + inventory + prepared splits ==="
rsync -avP results/hiv_env_mask/ "$BETTY:$REMOTE/results/hiv_env_mask/"
rsync -avP data/hiv_inventory.json "$BETTY:$REMOTE/data/"
rsync -avP data/hiv_temporal/ "$BETTY:$REMOTE/data/hiv_temporal/"
rsync -avP data/hiv_geo/ "$BETTY:$REMOTE/data/hiv_geo/"

if [[ "$SKIP_RAW" != "1" ]]; then
  echo "=== sync raw data/hiv (large) ==="
  rsync -avP data/hiv/ "$BETTY:$REMOTE/data/hiv/"
fi

if [[ "$SUBMIT" != "1" ]]; then
  echo "Sync only; SUBMIT=0"
  exit 0
fi

echo "=== submit on Betty ==="
ssh "$BETTY" bash -s <<'REMOTE_EOF'
set -euo pipefail
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints/hiv_temporal_v1 checkpoints/hiv_geo_v1 results/hiv_env_mask

# Ensure mask exists remotely
python scripts/build_hiv_env_v1v5_hotspot_mask.py \
  --out results/hiv_env_mask/mut_hotspot_mask_hiv_env_v1v5.pt

PIPE=$(sbatch --parsable scripts/slurm_hiv_pipeline.sh)
echo "PIPELINE_JOB=$PIPE"

TJOB=$(sbatch --parsable --dependency=afterok:$PIPE \
  --export=ALL,DATA=data/hiv_temporal,CKPT_DIR=checkpoints/hiv_temporal_v1,LAMBDA_BR=0 \
  --job-name=hiv_temp_v1 \
  scripts/slurm_hiv_train.sh)
echo "HIV_TEMPORAL_TRAIN=$TJOB"

GJOB=$(sbatch --parsable --dependency=afterok:$PIPE \
  --export=ALL,DATA=data/hiv_geo,CKPT_DIR=checkpoints/hiv_geo_v1,LAMBDA_BR=0 \
  --job-name=hiv_geo_v1 \
  scripts/slurm_hiv_train.sh)
echo "HIV_GEO_TRAIN=$GJOB"

echo "JOB_IDS pipeline=$PIPE temporal=$TJOB geo=$GJOB"
squeue -u nnori | head -20
REMOTE_EOF
