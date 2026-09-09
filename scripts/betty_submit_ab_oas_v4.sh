#!/bin/bash
# Antibody Recipe C only: cosine mut-head train → Rod.82.
# Does not overwrite paper ckpts, v2 SHM, v3 Thrifty, or ESM caches.
#
# On Betty: bash scripts/betty_submit_ab_oas_v4.sh

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs benchmarks/results/tables checkpoints/ab_oas_1m_v4_cosinehead

TRAIN_J=$(sbatch --parsable --qos=mig-max scripts/slurm_ab_oas_v4_cosinehead_train.sh)
echo "ab_oas_v4 cosinehead train -> $TRAIN_J"

ROD_J=$(sbatch --parsable --dependency=afterok:${TRAIN_J} \
  scripts/slurm_ab_oas_v4_cosinehead_rod82.sh)
echo "ab_oas_v4 rod82 -> $ROD_J"

python - <<PY
import json
from pathlib import Path
p = Path("benchmarks/results/tables/ab_oas_v4_cosinehead_job_ids.json")
p.write_text(json.dumps({
    "train": int("$TRAIN_J"),
    "rod82": int("$ROD_J"),
    "ckpt": "checkpoints/ab_oas_1m_v4_cosinehead",
    "write_as": "treesbm_ab_oas_v4",
    "mut_head": "cosine",
    "r0_backend": "esm2",
}, indent=2) + "\n")
print("Wrote", p)
PY

squeue -u "$USER" | head -30
