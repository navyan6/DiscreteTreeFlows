#!/bin/bash
#SBATCH --job-name=ab_oas_rod82
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/ab_oas_rod82_%j.log
#SBATCH --error=logs/ab_oas_rod82_%j.log
#
# After ab_oas_1m_v1 train: Rodriguez 82 TreeSBM rollout as treesbm_ab_oas
# (does not overwrite samples/treesbm or samples/treesbm_ab).
#
#   sbatch --dependency=afterok:<train_job> scripts/slurm_ab_oas_rod82_eval.sh
#
# If Rod.82 AA length > 160 (OAS ckpt), skip Rod and eval OAS test instead.

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

ROOT="${ROOT:-$HOME/DiscreteTreeFlows}"
cd "$ROOT"
PY="${PY:-/vast/home/n/nnori/.conda/envs/treesbm/bin/python}"
[ -x "$PY" ] || PY="${PY_FALLBACK:-/vast/projects/pranam/lab/nnori/.conda/envs/treesbm/bin/python}"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"

CKPT="${CKPT:-checkpoints/ab_oas_1m_v1/best.pt}"
CONFIG="${CONFIG:-antibody_benchmark/configs/full_ab_oas_treesbm.yaml}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-160}"
N_ROLLOUTS="${N_ROLLOUTS:-20}"
WRITE_AS="${WRITE_AS:-treesbm_ab_oas}"
SAMPLES_DIR="${SAMPLES_DIR:-antibody_benchmark/results/samples}"
OAS_RESULTS="${OAS_RESULTS:-antibody_benchmark/results_ab_oas}"

mkdir -p logs "$SAMPLES_DIR" "$OAS_RESULTS/samples" "$OAS_RESULTS/summary" \
  "$LABHOME/antibody_benchmark/logs"

echo "=== OAS TreeSBM Rodriguez 82 eval ==="
echo "host=$(hostname) date=$(date -Is) job=${SLURM_JOB_ID:-local}"
echo "ckpt=$CKPT write_as=$WRITE_AS max_seq_len=$MAX_SEQ_LEN"

if [[ ! -f "$CKPT" ]]; then
  echo "BLOCKER: missing $CKPT" >&2
  exit 1
fi

# Guard: never clobber existing Track C / Track B sample dirs
for keep in treesbm treesbm_ab; do
  if [[ "$WRITE_AS" == "$keep" ]]; then
    echo "REFUSING write-as=$WRITE_AS (would overwrite $keep)" >&2
    exit 2
  fi
done

rod_ok=$($PY - "$MAX_SEQ_LEN" <<'PY'
import json, sys
from pathlib import Path
L = int(sys.argv[1])
p = Path("antibody_benchmark/data/processed/benchmark_trees.jsonl")
if not p.is_file():
    print("no_jsonl")
    raise SystemExit(0)
n = 0
too = 0
mx = 0
for line in p.read_text().splitlines():
    if not line.strip():
        continue
    n += 1
    d = json.loads(line)
    aa = d.get("true_aa_sequences") or {}
    root = d.get("root_aa_sequence") or aa.get(d.get("root_id"), "") or ""
    lens = [len(root)] if root else []
    lens.extend(len(s) for s in aa.values() if s)
    if lens:
        mx = max(mx, max(lens))
        if max(lens) > L:
            too += 1
print(f"ok n={n} max_len={mx} n_over={too}")
PY
)
echo "rod_length_check: $rod_ok"

if [[ "$rod_ok" == no_jsonl ]] || [[ "$rod_ok" == *n_over=* && "$rod_ok" != *n_over=0* ]]; then
  echo "OAS ckpt max_seq_len=$MAX_SEQ_LEN cannot fairly eval Rodriguez 82 ($rod_ok)."
  echo "Falling back to OAS test clones (Table 6 Rod.82 OAS row stays empty)."
  mkdir -p results/ab_oas_1m_v1
  $PY scripts/eval_ab_maturation.py \
    --data data/ab_clones_1m/test \
    --methods plm_prior \
    --out results/ab_oas_1m_v1/eval_oas_test_fallback.json \
    --max-groups 170 \
    --K 100 || true
  echo "EVAL_TARGET=oas_test" > "$OAS_RESULTS/EVAL_TARGET.txt"
  echo "BLOCKER: Rod.82 length/alphabet mismatch vs OAS L=$MAX_SEQ_LEN. See $rod_ok"
  exit 0
fi

echo "EVAL_TARGET=rod82"
echo "rod82" > "$OAS_RESULTS/EVAL_TARGET.txt"

# Rollout into shared samples dir under a NEW model name
$PY antibody_benchmark/scripts/run_rollouts.py \
  --config "$CONFIG" \
  --models treesbm \
  --n-rollouts "$N_ROLLOUTS" \
  --out-dir "$SAMPLES_DIR" \
  --write-as "$WRITE_AS"

nsamp=$(find "$SAMPLES_DIR/$WRITE_AS" -name 'rollout_*.json' 2>/dev/null | wc -l | tr -d ' ')
echo "samples_written=$nsamp model=$WRITE_AS"
if [[ "$nsamp" -le 0 ]]; then
  echo "ERROR: no samples written under $SAMPLES_DIR/$WRITE_AS" >&2
  exit 1
fi

# OAS-only primary metrics (do not overwrite Track C summary/)
mkdir -p "$OAS_RESULTS/samples"
rm -f "$OAS_RESULTS/samples/$WRITE_AS"
ln -sfn "$(cd "$SAMPLES_DIR/$WRITE_AS" && pwd)" "$OAS_RESULTS/samples/$WRITE_AS"
$PY antibody_benchmark/scripts/evaluate.py --config "$CONFIG"

# Table 6 including OAS row (reads main samples dir; leaves treesbm / treesbm_ab intact)
$PY -u scripts/eval_ab_t6_from_samples.py \
  --trees antibody_benchmark/data/processed/benchmark_trees.jsonl \
  --samples-dir "$SAMPLES_DIR" \
  --models "${T6_MODELS:-thrifty,dasm_thrifty,cosine,treesbm,treesbm_ab,treesbm_ab_oas,${WRITE_AS}}" \
  --K 100 \
  --eps 1,2,3,5 \
  --out benchmarks/results/tables/table6_ab_track_c.json \
  --md-out benchmarks/results/tables/table6_ab_track_c.md

echo "Done: $(date -Is) samples=$nsamp target=rod82"
