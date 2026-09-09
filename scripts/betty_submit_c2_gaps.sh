#!/bin/bash
# Submit Appendix C.2 remaining cells on Betty (prefer login03; bash -lc for Slurm).
#
# Fills:
#   1) COVID TreeSBM Cov@500 — resume treesbm-only on K500 cache (missing NODE_0000179)
#   2) HIV Env geo pLM prior — coverage K500 (Cov@100/500, clade, min dist) + enrich (mut)
#   3) Ab pLM + treesbm_ab_oas — extend Rod.82 rollouts to 100 then Cov@500 + clade/min dist
#
# Usage (on Betty login03):
#   bash scripts/betty_submit_c2_gaps.sh
#   CANCEL_COVID=1 bash scripts/betty_submit_c2_gaps.sh   # cancel stuck/old 7628571 first

set -euo pipefail
cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results/tables checkpoints
export PATH=/cm/local/apps/slurm/current/bin:$PATH
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"

IDS=()
submit() {
  local name="$1"; shift
  local j
  j=$(sbatch --parsable --job-name="$name" "$@")
  echo "submitted $name -> $j"
  IDS+=("$name=$j")
}

OLD_COVID="${OLD_COVID_JOB:-7628571}"
if [[ "${CANCEL_COVID:-0}" == "1" ]]; then
  echo "Cancelling COVID Cov@500 job $OLD_COVID (resume treesbm-only next)…"
  scancel "$OLD_COVID" || true
  sleep 2
fi

# --- 1) COVID TreeSBM Cov@500 resume (cache has 4/5 roots) ---
submit "c2_covid_tsbm_k500" \
  --export=ALL,VIRUS=covid,METHODS="treesbm",K_MAX=500,K_STEP=50,MAX_ROOTS=5,N_LIST=16,CHECKPOINT=checkpoints/covid_v5_mutrec/best.pt,CACHE_DIR=benchmarks/results/coverage_cache_covid_N16_eabs_K500,OUT=benchmarks/results/coverage_curves_covid_N16_table5_eabs_K500.csv \
  --qos=mig-max --time=72:00:00 scripts/slurm_table5_coverage.sh

# --- 2) HIV geo pLM prior: coverage (all K + clade + min dist) ---
submit "c2_hivg_plm_cov" \
  --export=ALL,VIRUS=hiv_geo,METHODS="plm_prior",K_MAX=500,K_STEP=50,MAX_ROOTS=5,N_LIST=16,CACHE_DIR=benchmarks/results/coverage_cache_hiv_geo_N16_eabs_K500_plm,OUT=benchmarks/results/coverage_curves_hiv_geo_N16_table5_eabs_K500_plm.csv \
  --qos=mig-max --time=48:00:00 scripts/slurm_table5_coverage.sh

# Also K=100-only artifact path matching TreeSBM geo protocol (reuse if wanted)
submit "c2_hivg_plm_cov100" \
  --export=ALL,VIRUS=hiv_geo,METHODS="plm_prior",K_MAX=100,K_STEP=10,MAX_ROOTS=5,N_LIST=16,CACHE_DIR=benchmarks/results/coverage_cache_hiv_geo_N16_eabs_plm,OUT=benchmarks/results/coverage_curves_hiv_geo_N16_table5_eabs_plm.csv \
  --qos=mig-max --time=24:00:00 scripts/slurm_table5_coverage.sh

# HIV geo pLM enrich → Mut. recall (+ cons / aa|hit additives)
submit "c2_hivg_plm_enr" \
  --export=ALL,VIRUS=hiv_geo,METHODS="plm_prior",MAX_TREES=14,N=16,OUT=checkpoints/eval_table5_hiv_geo_plm_prior.json \
  --qos=mig-max --time=12:00:00 scripts/slurm_table5_baselines.sh

# --- 3) Ab: extend rollouts then eval Cov@500 + clade/min ---
submit "c2_ab_plm_roll" \
  --export=ALL,MODEL=plm_prior,N_ROLLOUTS=100 \
  --qos=mig-max --time=12:00:00 scripts/slurm_ab_c2_extend_rollouts.sh

submit "c2_ab_oas_roll" \
  --export=ALL,MODEL=treesbm,WRITE_AS=treesbm_ab_oas,CONFIG=antibody_benchmark/configs/full_ab_oas_treesbm.yaml,N_ROLLOUTS=100 \
  --qos=mig-max --time=12:00:00 scripts/slurm_ab_c2_extend_rollouts.sh

# Parse just-submitted roll job IDs for dependency
AB_PLM_J=""; AB_OAS_J=""
for kv in "${IDS[@]}"; do
  case "$kv" in
    c2_ab_plm_roll=*) AB_PLM_J="${kv#*=}" ;;
    c2_ab_oas_roll=*) AB_OAS_J="${kv#*=}" ;;
  esac
done
DEP="afterok:${AB_PLM_J}:${AB_OAS_J}"
# NOTE: do not put commas in --export values (Slurm splits on comma).
# slurm_ab_c2_eval.sh default MODELS=plm_prior,treesbm_ab_oas is correct.
submit "c2_ab_eval" \
  --dependency="$DEP" \
  --export=ALL,OUT=benchmarks/results/tables/table_c2_ab_metrics.json \
  scripts/slurm_ab_c2_eval.sh

printf '%s\n' "${IDS[@]}" | tee benchmarks/results/tables/table_c2_gap_job_ids.txt

python3 - <<'PY'
import json
from pathlib import Path
from datetime import datetime, timezone

lines = Path("benchmarks/results/tables/table_c2_gap_job_ids.txt").read_text().splitlines()
ids = {}
for ln in lines:
    if "=" in ln:
        k, v = ln.split("=", 1)
        ids[k] = v

# Merge into cov500 JSON (keep historical wave + new resumes)
cov_path = Path("benchmarks/results/tables/table_c2_cov500_job_ids.json")
prev = json.loads(cov_path.read_text()) if cov_path.is_file() else {}
payload = {
    **{k: v for k, v in prev.items() if k not in (
        "status_2026-08-16", "submitted", "fills_queued_2026-08-16", "jobs_gap_2026-08-16"
    )},
    "submitted": prev.get("submitted", "2026-08-15"),
    "gap_submit_2026-08-16": datetime.now(timezone.utc).isoformat(),
    "host": "betty_login03",
    "jobs_gap_2026-08-16": ids,
    "fills_queued_2026-08-16": {
        "c2_covid_tsbm_k500": "SARS-CoV-2 TreeSBM Cov@500 (resume NODE_0000179; cache 4/5)",
        "c2_hivg_plm_cov": "HIV Env geo pLM Cov@100+@500 + clade + min dist (K500 curve)",
        "c2_hivg_plm_cov100": "HIV Env geo pLM K=100 eabs CSV (protocol twin of TreeSBM geo)",
        "c2_hivg_plm_enr": "HIV Env geo pLM Mut. recall (enrich)",
        "c2_ab_plm_roll + c2_ab_oas_roll + c2_ab_eval": (
            "Ab pLM + treesbm_ab_oas Cov@500 + clade_recall + mean_min_edit "
            "(extend Rod.82 rollouts to 100, then CPU eval)"
        ),
    },
    "status_2026-08-16": {
        **(prev.get("status_2026-08-16") or {}),
        "covid_resume": f"queued treesbm-only → {ids.get('c2_covid_tsbm_k500')}",
        "hiv_plm": f"queued cov/enr → {ids.get('c2_hivg_plm_cov')}/{ids.get('c2_hivg_plm_enr')}",
        "ab_cov500_clade_mindist": (
            f"roll {ids.get('c2_ab_plm_roll')}/{ids.get('c2_ab_oas_roll')} "
            f"→ eval {ids.get('c2_ab_eval')}"
        ),
    },
}
# Preserve original virus→job map; annotate resume id
if "covid" in payload and ids.get("c2_covid_tsbm_k500"):
    payload["covid_resume_treesbm"] = ids["c2_covid_tsbm_k500"]
cov_path.write_text(json.dumps(payload, indent=2) + "\n")
print("wrote", cov_path)
print(json.dumps(ids, indent=2))
PY
