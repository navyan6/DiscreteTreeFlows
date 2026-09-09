#!/bin/bash
# Submit Table C.1 HIV baselines on Betty: NeutralBD / EmpiricalBD / AR / TreeSBM
# (topology via run_table) + Neutral/AR enrichment + Neutral/AR coverage.
#
# Usage (on Betty login):
#   bash scripts/betty_submit_hiv_c1_baselines.sh
#
# Prefer geo + temporal; N=16 first (24h wall). Uses existing
# checkpoints/hiv_{temporal,geo}_v1/best.pt and artreeformer_N16.nwk pools
# (topology-only adapted pools; BL fit on HIV train).

set -euo pipefail
cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results/tables

IDS=()

submit() {
  local name="$1"; shift
  local j
  j=$(sbatch --parsable --job-name="$name" "$@")
  echo "submitted $name -> $j"
  IDS+=("$name=$j")
}

# --- Topology / Tree-KL table (NeutralBD + EmpiricalBD + AR + TreeSBM) ---
# N=16 only to finish in 24h; sim_neutral for Tree-KL/Split-KL (empirical KL NaN).
for SPLIT in temporal geo; do
  DATA="data/hiv_${SPLIT}"
  CKPT="checkpoints/hiv_${SPLIT}_v1/best.pt"
  PARAMS="benchmarks/results/params_hiv_${SPLIT}.json"
  OUT="benchmarks/results/results_baselines_hiv_${SPLIT}.csv"
  TABLE_OUT="benchmarks/results/tables/hiv_${SPLIT}"
  submit "hiv_${SPLIT}_bl" --export=ALL,DATA="${DATA}/test",TRAIN_DATA="${DATA}/train",OUT="$OUT",TABLE_OUT="$TABLE_OUT",PARAMS="$PARAMS",N_LIST="16",K=20,M=20,MAX_ROOTS=14,MAX_SEQ_LEN=900,REGIMES="neutral",CKPT="$CKPT" \
    --qos=mig-max --time=24:00:00 scripts/slurm_baselines.sh "$CKPT"
done

# --- Sequence enrichment: NeutralBD + AR (V1–V5 antigenic) ---
for SPLIT in temporal geo; do
  submit "hiv_${SPLIT}_enr" --export=ALL,VIRUS="hiv_${SPLIT}",METHODS="neutral_bd artreeformer_adapted",MAX_TREES=14,N=16 \
    --qos=mig-max --time=12:00:00 scripts/slurm_table5_baselines.sh
done

# --- Coverage@100: Neutral + AR + TreeSBM (aa|hit / cons from curves) ---
for SPLIT in temporal geo; do
  submit "hiv_${SPLIT}_cov" --export=ALL,VIRUS="hiv_${SPLIT}",METHODS="neutral_bd artreeformer_adapted treesbm",N_LIST="16",MAX_ROOTS=5 \
    --qos=mig-max --time=12:00:00 scripts/slurm_table5_coverage.sh
done

printf '%s\n' "${IDS[@]}" | tee benchmarks/results/tables/hiv_c1_job_ids.txt
python3 - <<'PY'
import json
from pathlib import Path
from datetime import datetime, timezone
lines = Path("benchmarks/results/tables/hiv_c1_job_ids.txt").read_text().splitlines()
ids = {}
for ln in lines:
    if "=" in ln:
        k, v = ln.split("=", 1)
        ids[k] = v
Path("benchmarks/results/tables/hiv_c1_job_ids.json").write_text(json.dumps({
    "submitted_at": datetime.now(timezone.utc).isoformat(),
    "host": "betty",
    "jobs": ids,
    "note": "C.1 HIV: NeutralBD/EmpiricalBD/AR/TreeSBM baselines + Neutral/AR enrich + Neutral/AR/TreeSBM coverage. TreeSBM seq metrics already in eval_enrichment_hiv_*_mrs0.5_fullmetrics.json (7617937/39).",
}, indent=2) + "\n")
print("wrote hiv_c1_job_ids.json", ids)
PY
