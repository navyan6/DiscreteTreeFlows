#!/bin/bash
# Submit Appendix E.2 mut-head / stop-head train ablations + afterok enrich/coverage.
# 72h train wall (E.1 24h TIMEOUT lesson). --resume inside each ckpt dir.
set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints benchmarks/results/tables

submit_train() {
  local name="$1" ablate="$2"
  sbatch --parsable --qos=mig-max --time=72:00:00 --job-name="$name" \
    --export=ALL,ABLATE="$ablate",MODE=train \
    --output=logs/e2_ablate_%j.log --error=logs/e2_ablate_%j.log \
    scripts/slurm_e2_arch_ablations.sh
}

queue_eval() {
  local train_id="$1" ablate="$2" ckpt="$3" mode="$4" time="$5" jname="$6"
  sbatch --parsable --qos=mig-max --time="$time" --job-name="$jname" \
    --dependency="afterok:${train_id}" \
    --export=ALL,ABLATE="$ablate",MODE="$mode",CKPT="$ckpt" \
    --output=logs/e2_ablate_%j.log --error=logs/e2_ablate_%j.log \
    scripts/slurm_e2_arch_ablations.sh
}

J1=$(submit_train e2_no_mut no_mut_head)
J2=$(submit_train e2_no_stop no_stop_head)
echo "TRAIN: no_mut_head=$J1 no_stop_head=$J2"

E1=$(queue_eval "$J1" no_mut_head checkpoints/covid_e2_no_mut_head/best.pt enrich 12:00:00 e2en_nomut)
C1=$(queue_eval "$J1" no_mut_head checkpoints/covid_e2_no_mut_head/best.pt coverage 24:00:00 e2cv_nomut)
E2=$(queue_eval "$J2" no_stop_head checkpoints/covid_e2_no_stop_head/best.pt enrich 12:00:00 e2en_nostop)
C2=$(queue_eval "$J2" no_stop_head checkpoints/covid_e2_no_stop_head/best.pt coverage 24:00:00 e2cv_nostop)

echo "ENRICH: no_mut=$E1 no_stop=$E2"
echo "COV:    no_mut=$C1 no_stop=$C2"

python - <<PY
import json
from pathlib import Path
from datetime import datetime, timezone
ids = {
  "protocol": "covid_brazil_geo / covid_v5_mutrec recipe / mrs=0.5 / n_trees=20 / max_roots=5",
  "submitted": datetime.now(timezone.utc).isoformat(),
  "host": "betty_login03",
  "flags": {
    "no_mut_head": "--ablate-mut-head (train: log R_theta = log R0; lambda_mut=0)",
    "no_stop_head": "--ablate-stop-head (train: p_stop=0.5 constant; lambda_stop=0)",
  },
  "train": {
    "no_mut_head": int("$J1"),
    "no_stop_head": int("$J2"),
  },
  "enrich": {
    "no_mut_head": int("$E1"),
    "no_stop_head": int("$E2"),
  },
  "coverage": {
    "no_mut_head": int("$C1"),
    "no_stop_head": int("$C2"),
  },
  "ckpt_dirs": {
    "no_mut_head": "checkpoints/covid_e2_no_mut_head",
    "no_stop_head": "checkpoints/covid_e2_no_stop_head",
  },
  "timelimit": "72h (3-00:00:00) on train; --resume in each ckpt dir",
  "eval_metrics": "enrich: Mut/Cons/Ant/aa|hit/pLM NLL; coverage: Cov@100=coverage_obs_e2 @ K=100",
  "note": "Tree-KL not queued (saturates ln2). Stop head is L_stop-only at train; Alg. 4 does not sample p_stop.",
}
path = Path("benchmarks/results/tables/table_e2_job_ids.json")
path.write_text(json.dumps(ids, indent=2) + "\n")
print(json.dumps(ids, indent=2))
PY

squeue -u nnori | head -40
