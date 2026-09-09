#!/bin/bash
# Submit Appendix C.1 remaining cells on Betty (login03 preferred; bash -lc for Slurm PATH).
#
# Fills:
#   - Viral Neutral/AR pLM NLL (COVID / H3N2 / HIV geo)
#   - H3 TreeSBM pLM NLL (re-enrich; prior mrs0.5 JSON had null NLL)
#   - Antibody free-topo Tree-KL/Split-KL/RF/Quartet/W1 (slurm_baselines + OAS ckpt)
#   - Antibody Neutral/AR + TreeSBM pLM NLL
#
# Usage (on Betty login):
#   bash scripts/betty_submit_c1_gaps.sh

set -euo pipefail
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints benchmarks/results/tables
export PATH=/cm/local/apps/slurm/current/bin:$PATH

IDS=()
submit() {
  local name="$1"; shift
  local j
  j=$(sbatch --parsable --job-name="$name" "$@")
  echo "submitted $name -> $j"
  IDS+=("$name=$j")
}

# --- Viral Neutral / AR pLM NLL (free-topo gens + ESM-2 −PLL/pos) ---
submit "c1nll_covid" \
  --export=ALL,VIRUS=covid,METHODS="neutral_bd artreeformer_adapted",MAX_TREES=20,N=16,OUT=checkpoints/eval_table5_covid_baselines_plmnll.json \
  --qos=mig-max --time=12:00:00 scripts/slurm_table5_baselines.sh

submit "c1nll_h3" \
  --export=ALL,VIRUS=h3n2,METHODS="neutral_bd artreeformer_adapted",MAX_TREES=20,N=16,OUT=checkpoints/eval_table5_h3n2_baselines_plmnll.json \
  --qos=mig-max --time=12:00:00 scripts/slurm_table5_baselines.sh

submit "c1nll_hivg" \
  --export=ALL,VIRUS=hiv_geo,METHODS="neutral_bd artreeformer_adapted",MAX_TREES=14,N=16,OUT=checkpoints/eval_table5_hiv_geo_baselines_plmnll.json \
  --qos=mig-max --time=12:00:00 scripts/slurm_table5_baselines.sh

# --- H3 TreeSBM pLM NLL (prior enrich null) ---
submit "c1nll_h3ts" \
  --export=ALL,CKPT=checkpoints/h3n2_v3_lit_hotspot/best.pt,DATA=data/h3n2/test,MAX_SEQ_LEN=566,MRS=0.5,MAX_TREES=20,OUT=checkpoints/eval_enrichment_h3n2_v3_lit_hotspot_mrs0.5_plmnll.json,LIT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt \
  --qos=mig-max --time=12:00:00 scripts/slurm_c1_treesbm_plmnll.sh

# --- Antibody free-topo C.1 topology (NeutralBD / AR / TreeSBM OAS) ---
submit "c1ab_bl" \
  --export=ALL,DATA=data/ab_clones_1m/test,TRAIN_DATA=data/ab_clones_1m/train,OUT=benchmarks/results/results_baselines_ab_oas.csv,TABLE_OUT=benchmarks/results/tables/ab_oas,PARAMS=benchmarks/results/params_ab_oas.json,N_LIST="16",K=20,M=20,MAX_ROOTS=40,MAX_SEQ_LEN=160,REGIMES="neutral",CKPT=checkpoints/ab_oas_1m_v1/best.pt \
  --qos=mig-max --time=24:00:00 scripts/slurm_baselines.sh checkpoints/ab_oas_1m_v1/best.pt

# --- Antibody Neutral/AR pLM NLL ---
submit "c1ab_nll" \
  --export=ALL,VIRUS=ab_oas,METHODS="neutral_bd artreeformer_adapted",MAX_TREES=40,N=16,OUT=checkpoints/eval_table5_ab_oas_baselines_plmnll.json \
  --qos=mig-max --time=12:00:00 scripts/slurm_table5_baselines.sh

# --- Antibody TreeSBM (OAS) pLM NLL via free-gen enrich ---
submit "c1ab_tsnll" \
  --export=ALL,CKPT=checkpoints/ab_oas_1m_v1/best.pt,DATA=data/ab_clones_1m/test,MAX_SEQ_LEN=160,MRS=0.5,MAX_TREES=40,OUT=checkpoints/eval_enrichment_ab_oas_mrs0.5_plmnll.json,LIT_MASK=results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt \
  --qos=mig-max --time=12:00:00 scripts/slurm_c1_treesbm_plmnll.sh

printf '%s\n' "${IDS[@]}" | tee benchmarks/results/tables/table_c1_gap_job_ids.txt
python3 - <<'PY'
import json
from pathlib import Path
from datetime import datetime, timezone
lines = Path("benchmarks/results/tables/table_c1_gap_job_ids.txt").read_text().splitlines()
ids = {}
for ln in lines:
    if "=" in ln:
        k, v = ln.split("=", 1)
        ids[k] = v
payload = {
    "submitted_at": datetime.now(timezone.utc).isoformat(),
    "host": "betty",
    "protocol": "Appendix C.1 gap fill: Neutral/AR(+H3 TreeSBM) pLM NLL + Ab free-topo baselines",
    "jobs": ids,
    "fills": {
        "c1nll_covid": "SARS-CoV-2 Spike Neutral CTMC + Autoregressive tree-edit pLM NLL",
        "c1nll_h3": "Influenza HA H3N2 Neutral CTMC + Autoregressive tree-edit pLM NLL",
        "c1nll_hivg": "HIV Env Neutral CTMC + Autoregressive tree-edit pLM NLL (geo)",
        "c1nll_h3ts": "Influenza HA H3N2 TreeSBM pLM NLL",
        "c1ab_bl": "Antibody lineages Tree-KL/Split-KL/RF/Quartet/Branch W1 (all 3 methods)",
        "c1ab_nll": "Antibody lineages Neutral CTMC + Autoregressive tree-edit pLM NLL",
        "c1ab_tsnll": "Antibody lineages TreeSBM (OAS) pLM NLL",
    },
    "artifacts": {
        "c1nll_covid": "checkpoints/eval_table5_covid_baselines_plmnll.json",
        "c1nll_h3": "checkpoints/eval_table5_h3n2_baselines_plmnll.json",
        "c1nll_hivg": "checkpoints/eval_table5_hiv_geo_baselines_plmnll.json",
        "c1nll_h3ts": "checkpoints/eval_enrichment_h3n2_v3_lit_hotspot_mrs0.5_plmnll.json",
        "c1ab_bl": "benchmarks/results/tables/ab_oas/table_empirical_N16.csv (+ sim_neutral)",
        "c1ab_nll": "checkpoints/eval_table5_ab_oas_baselines_plmnll.json",
        "c1ab_tsnll": "checkpoints/eval_enrichment_ab_oas_mrs0.5_plmnll.json",
    },
}
Path("benchmarks/results/tables/table_c1_job_ids.json").write_text(
    json.dumps(payload, indent=2) + "\n"
)
print("wrote table_c1_job_ids.json", ids)
PY
