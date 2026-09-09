#!/bin/bash
# Antibody Recipe B only: Thrifty Q0 precompute → train → Rod.82.
# Does not overwrite paper ckpts, v2 SHM, or ESM group_*_ref_rates.pt.
#
# On Betty: bash scripts/betty_submit_ab_oas_v3.sh

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs benchmarks/results/tables checkpoints/ab_oas_1m_v3_thrifty

PRE_J=$(sbatch --parsable --qos=mig-max scripts/slurm_ab_oas_v3_thrifty_precompute.sh)
echo "ab_oas_v3 thrifty precompute -> $PRE_J"

TRAIN_J=$(sbatch --parsable --qos=mig-max --dependency=afterok:${PRE_J} \
  scripts/slurm_ab_oas_v3_thrifty_train.sh)
echo "ab_oas_v3 thrifty train -> $TRAIN_J"

ROD_J=$(sbatch --parsable --dependency=afterok:${TRAIN_J} \
  scripts/slurm_ab_oas_v3_thrifty_rod82.sh)
echo "ab_oas_v3 rod82 -> $ROD_J"

python - <<PY
import json
from pathlib import Path
p = Path("benchmarks/results/tables/ab_oas_v3_thrifty_job_ids.json")
p.write_text(json.dumps({
    "precompute": int("$PRE_J"),
    "train": int("$TRAIN_J"),
    "rod82": int("$ROD_J"),
    "ckpt": "checkpoints/ab_oas_1m_v3_thrifty",
    "write_as": "treesbm_ab_oas_v3",
    "r0_backend": "thrifty_aa",
}, indent=2) + "\n")
print("Wrote", p)
PY

squeue -u "$USER" | head -30
