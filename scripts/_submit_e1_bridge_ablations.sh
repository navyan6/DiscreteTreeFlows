#!/bin/bash
# Submit Appendix E.1 bridge-matching train ablations + afterok enrich/coverage.
set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints benchmarks/results/tables

submit_train() {
  local name="$1" ablate="$2"
  sbatch --qos=mig-max --time=24:00:00 --job-name="$name" \
    --export=ALL,ABLATE="$ablate",MODE=train \
    --output=logs/e1_ablate_%j.log --error=logs/e1_ablate_%j.log \
    scripts/slurm_e1_bridge_ablations.sh | awk '{print $4}'
}

queue_eval() {
  local train_id="$1" ablate="$2" ckpt="$3" mode="$4" time="$5" jname="$6"
  sbatch --qos=mig-max --time="$time" --job-name="$jname" \
    --dependency="afterok:${train_id}" \
    --export=ALL,ABLATE="$ablate",MODE="$mode",CKPT="$ckpt" \
    --output=logs/e1_ablate_%j.log --error=logs/e1_ablate_%j.log \
    scripts/slurm_e1_bridge_ablations.sh | awk '{print $4}'
}

J1=$(submit_train e1_term_only terminal_only)
J2=$(submit_train e1_no_doob no_doob)
J3=$(submit_train e1_no_term no_terminal)
echo "TRAIN: terminal_only=$J1 no_doob=$J2 no_terminal=$J3"

E1=$(queue_eval "$J1" terminal_only checkpoints/covid_e1_terminal_only/best.pt enrich 12:00:00 e1en_term)
C1=$(queue_eval "$J1" terminal_only checkpoints/covid_e1_terminal_only/best.pt coverage 24:00:00 e1cv_term)
E2=$(queue_eval "$J2" no_doob checkpoints/covid_e1_no_doob/best.pt enrich 12:00:00 e1en_nodoob)
C2=$(queue_eval "$J2" no_doob checkpoints/covid_e1_no_doob/best.pt coverage 24:00:00 e1cv_nodoob)
E3=$(queue_eval "$J3" no_terminal checkpoints/covid_e1_no_terminal/best.pt enrich 12:00:00 e1en_noterm)
C3=$(queue_eval "$J3" no_terminal checkpoints/covid_e1_no_terminal/best.pt coverage 24:00:00 e1cv_noterm)

echo "ENRICH: term=$E1 no_doob=$E2 no_term=$E3"
echo "COV:    term=$C1 no_doob=$C2 no_term=$C3"

python - <<PY
import json
from pathlib import Path
ids = {
  "protocol": "covid_brazil_geo / covid_v5_mutrec recipe / mrs=0.5 / n_trees=20 / max_roots=5",
  "submitted": "2026-08-15",
  "flags": {
    "terminal_only": "--ablate-terminal-only (train: force t->1 CE target)",
    "no_doob": "--ablate-doob (train: KL(softmax(R0)||R_theta), no Doob h-transform)",
    "no_terminal": "--ablate-terminal-consistency (train: lambda_cons=0)",
    "reference_only": "T8 --ablate-bridge on covid_v5_mutrec (already done)",
    "full": "T8 full on covid_v5_mutrec (already done)",
  },
  "train": {
    "terminal_only": int("$J1"),
    "no_doob": int("$J2"),
    "no_terminal": int("$J3"),
  },
  "enrich": {
    "terminal_only": int("$E1"),
    "no_doob": int("$E2"),
    "no_terminal": int("$E3"),
  },
  "coverage": {
    "terminal_only": int("$C1"),
    "no_doob": int("$C2"),
    "no_terminal": int("$C3"),
  },
  "ckpt_dirs": {
    "terminal_only": "checkpoints/covid_e1_terminal_only",
    "no_doob": "checkpoints/covid_e1_no_doob",
    "no_terminal": "checkpoints/covid_e1_no_terminal",
  },
  "note": "prior lambda_br=0 retrain 7539981 TIMEOUT 24h; these may also leave best.pt mid-run; afterok enrich/cov depend on train success",
}
path = Path("benchmarks/results/tables/table_e1_job_ids.json")
path.write_text(json.dumps(ids, indent=2) + "\n")
print(json.dumps(ids, indent=2))
PY

squeue -u nnori | head -40
