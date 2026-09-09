#!/bin/bash
# Appendix E.3 — sampling hyperparameter sensitivity (COVID Brazil / covid_v5_mutrec).
#
# One-factor-at-a-time (not 3^3 factorial). Defaults: K_max=100, T=1.0,
# max_leaves = production N-adapter (≈ N*1.6+2; documented as ~400 gen default).
#
# Grid values:
#   K ∈ {50, 100, 500}           — from existing T5 eabs + C.2 K500 (7635911)
#   temperature ∈ {0.5, 1.0, 1.5} — submit 0.5 & 1.5; 1.0 = Table5 treesbm baseline
#   max_leaves ∈ {100, 400, 800}  — submit 100 & 800; 400 ≈ production default row
#
# Usage (Betty login03):
#   bash scripts/betty_submit_e3.sh

set -euo pipefail
cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results/tables
export PATH=/cm/local/apps/slurm/current/bin:$PATH
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"

CKPT=checkpoints/covid_v5_mutrec/best.pt
IDS=()
submit() {
  local name="$1"; shift
  local j
  j=$(sbatch --parsable --job-name="$name" "$@")
  echo "submitted $name -> $j"
  IDS+=("$name=$j")
}

for T in 0.5 1.5; do
  tag=$(echo "$T" | tr '.' 'p')
  submit "e3_temp_${tag}" \
    --export=ALL,VIRUS=covid,METHODS="treesbm",K_MAX=100,K_STEP=50,MAX_ROOTS=5,N_LIST=16,CHECKPOINT="$CKPT",SITE_TEMPERATURE="$T",CACHE_DIR="benchmarks/results/coverage_cache_covid_e3_temp_${tag}",OUT="benchmarks/results/coverage_curves_covid_N16_e3_temp_${tag}.csv" \
    --qos=mig-max --time=48:00:00 scripts/slurm_table5_coverage.sh
done

for ML in 100 800; do
  submit "e3_maxl_${ML}" \
    --export=ALL,VIRUS=covid,METHODS="treesbm",K_MAX=100,K_STEP=50,MAX_ROOTS=5,N_LIST=16,CHECKPOINT="$CKPT",SITE_TEMPERATURE=1.0,MAX_LEAVES="$ML",CACHE_DIR="benchmarks/results/coverage_cache_covid_e3_maxl_${ML}",OUT="benchmarks/results/coverage_curves_covid_N16_e3_maxl_${ML}.csv" \
    --qos=mig-max --time=48:00:00 scripts/slurm_table5_coverage.sh
done

printf '%s\n' "${IDS[@]}" | tee benchmarks/results/tables/table_e3_gap_job_ids.txt

python3 - <<'PY'
import json
from pathlib import Path
from datetime import datetime, timezone

lines = Path("benchmarks/results/tables/table_e3_gap_job_ids.txt").read_text().splitlines()
ids = {}
for ln in lines:
    if "=" in ln:
        k, v = ln.split("=", 1)
        ids[k] = v

payload = {
    "submitted": datetime.now(timezone.utc).isoformat(),
    "host": "betty_login03",
    "protocol": {
        "virus": "covid",
        "ckpt": "checkpoints/covid_v5_mutrec/best.pt",
        "roots": "Table5 Brazil groups 4/5/7/9/10",
        "N": 16,
        "primary_metric": "coverage_obs_e2 @ K (same as Table5 / C.2)",
        "design": "one-factor-at-a-time (not full factorial)",
        "K": [50, 100, 500],
        "temperature": [0.5, 1.0, 1.5],
        "max_leaves": [100, 400, 800],
        "defaults_held": {
            "K_max": 100,
            "site_temperature": 1.0,
            "max_leaves": "adapter N*cushion+2 (production) / documented 400 default",
        },
        "mapped_not_resubmitted": {
            "K_50_100_T1": "coverage_curves_covid_N16_table5_eabs.csv treesbm",
            "K_500": "C.2 job 7635911 → coverage_curves_covid_N16_table5_eabs_K500.csv",
            "T_1.0": "same as Table5 treesbm row",
            "max_leaves_400": "production generate_k default; adapter path ≈ Table5",
        },
    },
    "jobs_submitted": ids,
    "artifacts": {
        "K_baseline": "benchmarks/results/coverage_curves_covid_N16_table5_eabs.csv",
        "K500": "benchmarks/results/coverage_curves_covid_N16_table5_eabs_K500.csv",
        "temp_0p5": "benchmarks/results/coverage_curves_covid_N16_e3_temp_0p5.csv",
        "temp_1p5": "benchmarks/results/coverage_curves_covid_N16_e3_temp_1p5.csv",
        "maxl_100": "benchmarks/results/coverage_curves_covid_N16_e3_maxl_100.csv",
        "maxl_800": "benchmarks/results/coverage_curves_covid_N16_e3_maxl_800.csv",
    },
}
Path("benchmarks/results/tables/table_e3_job_ids.json").write_text(
    json.dumps(payload, indent=2) + "\n"
)
print(json.dumps(ids, indent=2))
PY
