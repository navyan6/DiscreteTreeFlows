#!/bin/bash
# Submit Rod.82 Neutral SHM / pLM / AR rollouts + eval for paper Table 5 Ab.
# Usage (on Betty): bash antibody_benchmark/scripts/betty_submit_t5_baselines.sh

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
ROOT="${ROOT:-$HOME/DiscreteTreeFlows}"
cd "$ROOT"
mkdir -p logs antibody_benchmark/results/samples antibody_benchmark/results/smoke_reports

# Neutral SHM — CPU, JC69
J_NEUT=$(MODEL=neutral_shm sbatch --parsable \
  --job-name=ab_neutral_shm \
  --partition=genoa-std-mem \
  --cpus-per-task=4 \
  --mem-per-cpu=5632M \
  --time=6:00:00 \
  --output=logs/ab_neutral_shm_%j.log \
  --error=logs/ab_neutral_shm_%j.log \
  antibody_benchmark/scripts/slurm_ab_full_rollout.sh)
echo "SUBMITTED neutral_shm -> $J_NEUT"

# AR tree-edit — CPU (pyvolve + pool prune)
J_AR=$(MODEL=ar_tree_edit sbatch --parsable \
  --job-name=ab_ar_tree_edit \
  --partition=genoa-std-mem \
  --cpus-per-task=4 \
  --mem-per-cpu=5632M \
  --time=6:00:00 \
  --output=logs/ab_ar_tree_edit_%j.log \
  --error=logs/ab_ar_tree_edit_%j.log \
  antibody_benchmark/scripts/slurm_ab_full_rollout.sh)
echo "SUBMITTED ar_tree_edit -> $J_AR"

# pLM prior — GPU ESM
J_PLM=$(MODEL=plm_prior sbatch --parsable \
  --job-name=ab_plm_prior \
  --partition=b200-mig45 \
  --qos=mig-max \
  --gres=gpu:1 \
  --cpus-per-task=6 \
  --time=12:00:00 \
  --output=logs/ab_plm_prior_%j.log \
  --error=logs/ab_plm_prior_%j.log \
  antibody_benchmark/scripts/slurm_ab_full_rollout.sh)
echo "SUBMITTED plm_prior -> $J_PLM"

# Eval after all three (also recomputes existing models)
J_EVAL=$(sbatch --parsable \
  --job-name=ab_t5_eval_metrics \
  --partition=genoa-std-mem \
  --cpus-per-task=4 \
  --mem-per-cpu=5632M \
  --time=2:00:00 \
  --dependency=afterok:${J_NEUT}:${J_AR}:${J_PLM} \
  --output=logs/ab_t5_eval_metrics_%j.log \
  --error=logs/ab_t5_eval_metrics_%j.log \
  --wrap="cd ~/DiscreteTreeFlows && export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:\$PATH && export PYTHONPATH=~/DiscreteTreeFlows && python scripts/eval_ab_t6_from_samples.py --samples-dir antibody_benchmark/results/samples --models neutral_shm,plm_prior,ar_tree_edit,cosine,treesbm,treesbm_ab,treesbm_ab_oas,thrifty,dasm_thrifty --out benchmarks/results/tables/table6_ab_track_c.json --md-out benchmarks/results/tables/table6_ab_track_c.md")
echo "SUBMITTED eval afterok -> $J_EVAL (dep ${J_NEUT}:${J_AR}:${J_PLM})"

python3 - <<PY
import json
ids = {
  "neutral_shm": "$J_NEUT",
  "ar_tree_edit": "$J_AR",
  "plm_prior": "$J_PLM",
  "eval": "$J_EVAL",
}
path = "antibody_benchmark/results/smoke_reports/t5_baseline_submit_ids.json"
with open(path, "w") as f:
    json.dump(ids, f, indent=2)
    f.write("\n")
print(json.dumps(ids, indent=2))
print("wrote", path)
PY
