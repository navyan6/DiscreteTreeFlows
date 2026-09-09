#!/usr/bin/env bash
# Submit full 82×N=20 antibody benchmark jobs on Betty (run ON Betty after sync).
set -euo pipefail
# Prefer module env (Betty login nodes need this for Slurm DNS/config)
if type module >/dev/null 2>&1; then
  module load slurm 2>/dev/null || true
fi
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
ROOT="${ROOT:-$HOME/DiscreteTreeFlows}"
cd "$ROOT"
mkdir -p "$LABHOME/antibody_benchmark/logs"

SBATCH_SH=antibody_benchmark/scripts/slurm_ab_full_rollout.sh
chmod +x "$SBATCH_SH" antibody_benchmark/scripts/betty_submit_full.sh 2>/dev/null || true

# Preflight
test -f antibody_benchmark/data/processed/benchmark_trees.jsonl
n_trees=$(wc -l < antibody_benchmark/data/processed/benchmark_trees.jsonl | tr -d ' ')
echo "n_trees=$n_trees (expect 82)"
test -f checkpoints/best.pt
CKPT_LAB="$LABHOME/antibody_benchmark/cosine_ckpts/cosine_dasm.ckpt"
test -f "$CKPT_LAB" || test -f antibody_benchmark/data/raw/cosine/checkpoints/cosine_dasm.ckpt

declare -A JOBS=()

# CPU jobs: Thrifty, DASM+Thrifty (24h); TreeSBM (48h — ESM-heavy)
for MODEL in thrifty dasm_thrifty; do
  jid=$(MODEL="$MODEL" sbatch --parsable \
    --job-name="ab_${MODEL}" \
    --time=24:00:00 \
    "$SBATCH_SH")
  JOBS[$MODEL]=$jid
  echo "SUBMITTED $MODEL -> $jid (CPU genoa-std-mem, 24h)"
done

jid=$(MODEL=treesbm sbatch --parsable \
  --job-name=ab_treesbm \
  --time=48:00:00 \
  "$SBATCH_SH")
JOBS[treesbm]=$jid
echo "SUBMITTED treesbm -> $jid (CPU genoa-std-mem, 48h)"

# CoSiNE GPU (same partition/qos as successful smoke); 48h wall for 82×20
jid=$(MODEL=cosine sbatch --parsable \
  --job-name=ab_cosine \
  --partition=b200-mig45 \
  --qos=mig-max \
  --gres=gpu:1 \
  --cpus-per-task=6 \
  --mem-per-cpu=8G \
  --time=48:00:00 \
  "$SBATCH_SH")
JOBS[cosine]=$jid
echo "SUBMITTED cosine -> $jid (GPU b200-mig45 mig-max, 48h)"

# Evaluate after all rollouts (dependency)
dep="${JOBS[thrifty]}:${JOBS[dasm_thrifty]}:${JOBS[treesbm]}:${JOBS[cosine]}"
EVAL_SH="$LABHOME/antibody_benchmark/ab_full_evaluate.sbatch"
cat > "$EVAL_SH" <<EOF
#!/bin/bash
#SBATCH -J ab_eval
#SBATCH -p genoa-std-mem
#SBATCH -c 4
#SBATCH --mem-per-cpu=5632M
#SBATCH -t 04:00:00
#SBATCH -o $LABHOME/antibody_benchmark/logs/ab_eval_%j.out
#SBATCH -e $LABHOME/antibody_benchmark/logs/ab_eval_%j.err
set -euo pipefail
cd "\$HOME/DiscreteTreeFlows"
PY=/vast/projects/pranam/lab/nnori/.conda/envs/treesbm/bin/python
export PYTHONPATH="\$PWD\${PYTHONPATH:+:\$PYTHONPATH}"
"\$PY" antibody_benchmark/scripts/evaluate.py --config antibody_benchmark/configs/full.yaml
echo Done: \$(date -Is)
EOF
ejid=$(sbatch --parsable --dependency=afterok:${dep} "$EVAL_SH")
JOBS[evaluate]=$ejid
echo "SUBMITTED evaluate -> $ejid (afterok:$dep)"

OUT="$LABHOME/antibody_benchmark/full_submit_ids.json"
python3 - <<PY
import json, os
jobs = {
  "thrifty": "${JOBS[thrifty]}",
  "dasm_thrifty": "${JOBS[dasm_thrifty]}",
  "treesbm": "${JOBS[treesbm]}",
  "cosine": "${JOBS[cosine]}",
  "evaluate": "${JOBS[evaluate]}",
}
meta = {
  "submitted_utc": __import__("datetime").datetime.utcnow().isoformat() + "Z",
  "n_trees": int("$n_trees"),
  "n_rollouts": 20,
  "config": "antibody_benchmark/configs/full.yaml",
  "samples_dir": "antibody_benchmark/results/samples/<model>/<family>/rollout_*.json",
  "summary_dir": "antibody_benchmark/results/summary/",
  "logs_dir": "$LABHOME/antibody_benchmark/logs/",
  "jobs": jobs,
}
path = "$OUT"
open(path, "w").write(json.dumps(meta, indent=2))
print(json.dumps(meta, indent=2))
PY

squeue -u nnori | head -20
echo "WROTE $OUT"
