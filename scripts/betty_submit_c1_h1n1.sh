#!/bin/bash
# Submit Appendix C.1 Influenza HA H1N1 cells on Betty (login03 preferred).
#
# Replaces HIV Env in Table C.1. Best TreeSBM = geo-split Table 5 ckpt
# checkpoints/h1n1_v2_lit_hotspot/best.pt (NOT leafholdout / temporal).
#
# Fills:
#   - Topology N=16 empirical (NeutralBD / AR / TreeSBM): RF / Quartet / W1
#     + sim_neutral Tree-KL / Split-KL
#   - Neutral/AR pLM NLL
#   - TreeSBM pLM NLL (prior mrs0.5 enrich JSON had no plm_nll)
#
# Analogous to C.1 NLL jobs 7635749–52 + Ab/HIV free-topo baselines.
#
# Usage (on Betty login):
#   bash scripts/betty_submit_c1_h1n1.sh

set -euo pipefail
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints benchmarks/results/tables/h1n1
export PATH=/cm/local/apps/slurm/current/bin:$PATH

IDS=()
submit() {
  local name="$1"; shift
  local j
  j=$(sbatch --parsable --job-name="$name" "$@")
  echo "submitted $name -> $j"
  IDS+=("$name=$j")
}

# --- Topology / Tree-KL table (NeutralBD + AR + TreeSBM geo lit-hotspot) ---
# N=16 only; sim_neutral for Tree-KL/Split-KL (empirical KL NaN).
# MAX_ROOTS=40 matches geo test tree count (table2_datasets.md).
submit "c1h1_bl" \
  --export=ALL,DATA=data/h1n1/test,TRAIN_DATA=data/h1n1/train,OUT=benchmarks/results/results_baselines_h1n1.csv,TABLE_OUT=benchmarks/results/tables/h1n1,PARAMS=benchmarks/results/params_h1n1.json,N_LIST="16",K=20,M=20,MAX_ROOTS=40,MAX_SEQ_LEN=566,REGIMES="neutral",CKPT=checkpoints/h1n1_v2_lit_hotspot/best.pt \
  --qos=mig-max --time=24:00:00 scripts/slurm_baselines.sh checkpoints/h1n1_v2_lit_hotspot/best.pt

# --- Neutral / AR pLM NLL (free-topo gens + ESM-2 −PLL/pos) ---
submit "c1nll_h1" \
  --export=ALL,VIRUS=h1n1,METHODS="neutral_bd artreeformer_adapted",MAX_TREES=20,N=16,OUT=checkpoints/eval_table5_h1n1_baselines_plmnll.json \
  --qos=mig-max --time=12:00:00 scripts/slurm_table5_baselines.sh

# --- TreeSBM pLM NLL (prior mrs0.5 enrich had no plm_nll) ---
submit "c1nll_h1ts" \
  --export=ALL,CKPT=checkpoints/h1n1_v2_lit_hotspot/best.pt,DATA=data/h1n1/test,MAX_SEQ_LEN=566,MRS=0.5,MAX_TREES=20,OUT=checkpoints/eval_enrichment_h1n1_v2_lit_hotspot_mrs0.5_plmnll.json,LIT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_h1_nmicrobiol_lit.pt \
  --qos=mig-max --time=12:00:00 scripts/slurm_c1_treesbm_plmnll.sh

printf '%s\n' "${IDS[@]}" | tee benchmarks/results/tables/table_c1_h1n1_gap_job_ids.txt
python3 - <<'PY'
import json
from pathlib import Path
from datetime import datetime, timezone
lines = Path("benchmarks/results/tables/table_c1_h1n1_gap_job_ids.txt").read_text().splitlines()
ids = {}
for ln in lines:
    if "=" in ln:
        k, v = ln.split("=", 1)
        ids[k] = v
payload = {
    "submitted_at": datetime.now(timezone.utc).isoformat(),
    "host": "betty",
    "protocol": "Appendix C.1 H1N1 replace HIV Env: topo N=16 + Neutral/AR/TreeSBM pLM NLL",
    "best_treesbm": {
        "ckpt": "checkpoints/h1n1_v2_lit_hotspot/best.pt",
        "split": "geo (data/h1n1) — same as Table 5 / C.3; NOT h1n1_leafholdout or h1n1_temporal",
        "why": "Table 5 viral forecasting + C.3 horizon locked to this ckpt (ep74). leafholdout/temporal unused for paper T5.",
    },
    "jobs": ids,
    "fills": {
        "c1h1_bl": "Influenza HA H1N1 Tree-KL/Split-KL/RF/Quartet/Branch W1 (Neutral CTMC / AR / TreeSBM)",
        "c1nll_h1": "Influenza HA H1N1 Neutral CTMC + Autoregressive tree-edit pLM NLL",
        "c1nll_h1ts": "Influenza HA H1N1 TreeSBM pLM NLL",
    },
    "artifacts": {
        "c1h1_bl": "benchmarks/results/results_baselines_h1n1.csv + tables/h1n1/table_empirical_N16.csv (+ sim_neutral)",
        "c1nll_h1": "checkpoints/eval_table5_h1n1_baselines_plmnll.json",
        "c1nll_h1ts": "checkpoints/eval_enrichment_h1n1_v2_lit_hotspot_mrs0.5_plmnll.json",
    },
    "analogous_to": {
        "c1nll_h1": "7635749 covid / 7635750 h3 / 7635751 hiv_geo",
        "c1nll_h1ts": "7635752 h3 TreeSBM",
        "c1h1_bl": "7635753 Ab free-topo / 7626610 HIV geo topo",
    },
    "prior_artifacts_no_nll": {
        "eval_table5_h1n1_baselines.json": "job 7539950 — plm_prior + AR only; no plm_nll field",
        "eval_enrichment_h1n1_v2_lit_hotspot_mrs0.5_fullmetrics.json": "no plm_nll",
        "results_baselines_h1n1.csv": "absent (no prior topo suite)",
    },
}
Path("benchmarks/results/tables/table_c1_h1n1_job_ids.json").write_text(
    json.dumps(payload, indent=2) + "\n"
)
print("wrote table_c1_h1n1_job_ids.json", ids)
PY
