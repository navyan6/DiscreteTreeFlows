#!/bin/bash
# 1) Fixed-topo eval on paper ckpts (COVID/H1/H3/HIV geo)
# 2) COVID codon-GY94 Q0 train → afterok fixed-topo eval of new ckpt
# Does not overwrite paper / mutlin / OAS ckpts.
#
# On Betty: bash scripts/betty_submit_fixed_topo_codon.sh

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs checkpoints/viral_eval_sweep benchmarks/results/tables

ids=()
for v in covid h1n1 h3n2 hiv_geo; do
  J=$(sbatch --parsable --qos=mig-max --export=ALL,VIRUS="$v" scripts/slurm_eval_fixed_topo.sh)
  echo "fixed_topo $v -> $J"
  ids+=("$J")
done

TRAIN_J=$(sbatch --parsable --qos=mig-max scripts/slurm_covid_codon_q0_train.sh)
echo "covid_v8_codon_q0 train -> $TRAIN_J"

EVAL8=$(sbatch --parsable --qos=mig-max --dependency=afterok:${TRAIN_J} \
  --export=ALL,VIRUS=covid,CKPT=checkpoints/covid_v8_codon_q0/best.pt,OUT=checkpoints/viral_eval_sweep/covid_v8_codon_q0_fixed_topo.json \
  scripts/slurm_eval_fixed_topo.sh)
echo "covid_v8 fixed_topo -> $EVAL8"

python - <<PY
import json
from pathlib import Path
p = Path("benchmarks/results/tables/fixed_topo_codon_job_ids.json")
p.write_text(json.dumps({
    "fixed_topo_paper": {
        "covid": int("${ids[0]}"),
        "h1n1": int("${ids[1]}"),
        "h3n2": int("${ids[2]}"),
        "hiv_geo": int("${ids[3]}"),
    },
    "codon_q0_train": int("$TRAIN_J"),
    "codon_q0_eval": int("$EVAL8"),
    "ckpt_new": "checkpoints/covid_v8_codon_q0",
    "locked_untouched": [
        "covid_v5_mutrec", "h1n1_v2_lit_hotspot", "h3n2_v3_lit_hotspot",
        "hiv_geo_v1", "ab_oas_1m_v1",
    ],
}, indent=2) + "\n")
print("Wrote", p)
PY
squeue -u "$USER" | head -20
