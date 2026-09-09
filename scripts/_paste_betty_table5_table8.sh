#!/bin/bash
# Paste / run from laptop: sync Table 5/8 scripts to Betty and submit jobs.
# No git push. Requires working SSH (kinit + Duo if needed).
set -euo pipefail
HOST=nnori@login.betty.parcc.upenn.edu
REMOTE=~/DiscreteTreeFlows

rsync -avz \
  scripts/eval_table5_baselines.py \
  scripts/eval_single_tree.py \
  scripts/eval_evescape_enrichment.py \
  scripts/slurm_table5_baselines.sh \
  scripts/slurm_table5_coverage.sh \
  scripts/slurm_table8_ablations.sh \
  scripts/slurm_covid_train_mut_recovery.sh \
  benchmarks/coverage_curves.py \
  benchmarks/results/tables/table5_viral_forecasting.csv \
  benchmarks/results/tables/table5_viral_forecasting.md \
  benchmarks/results/tables/table8_ablations.md \
  "$HOST:$REMOTE/"

# Flatten into correct remote paths
ssh -o BatchMode=yes "$HOST" 'bash -s' <<'REMOTE'
set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:/opt/slurm/bin:$PATH
export SLURM_CONF=/cm/shared/apps/slurm/etc/slurm/slurm.conf
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints benchmarks/results/tables scripts

# rsync may have dropped files in home if relative; ensure paths
# (rsync with trailing paths above puts files preserving structure under REMOTE)

echo "==== inventory ===="
ls -la scripts/slurm_table5_*.sh scripts/slurm_table8_*.sh scripts/eval_table5_baselines.py 2>&1 | head
ls benchmarks/external_pools/sampled/ 2>&1 | head
test -f checkpoints/covid_v5_mutrec/best.pt && echo covid_v5_ok || echo covid_v5_MISSING
test -f checkpoints/covid_v7_pmc_hotspot/best.pt && echo covid_v7_ok || echo covid_v7_MISSING
test -f checkpoints/h1n1_v2_lit_hotspot/best.pt && echo h1_ok || echo h1_MISSING
test -f benchmarks/external_pools/sampled/artreeformer_N16.nwk && echo ar_pool_ok || echo ar_pool_MISSING
test -d data/covid/test && ls data/covid/test/*rooted.nwk 2>/dev/null | wc -l || echo no_covid_test
test -d data/h1n1/test && ls data/h1n1/test/*rooted.nwk 2>/dev/null | wc -l || echo no_h1_test

JOBS=()
submit() {
  local j
  j=$(sbatch --qos=mig-max "$@" | awk '{print $4}')
  echo "SUBMITTED $j :: $*"
  JOBS+=("$j")
}

# Table 5 coverage (COVID + H1)
submit --job-name=t5_cov_covid --export=ALL,VIRUS=covid scripts/slurm_table5_coverage.sh
submit --job-name=t5_cov_h1 --export=ALL,VIRUS=h1n1 scripts/slurm_table5_coverage.sh

# Table 5 baselines — prefer PLM+AR; if AR pool missing, PLM-only still runs
if [[ -f benchmarks/external_pools/sampled/artreeformer_N16.nwk ]]; then
  METH="plm_prior artreeformer_adapted"
else
  METH="plm_prior"
  echo "WARN: artreeformer_N16.nwk missing — submitting PLM-only baselines"
fi
submit --job-name=t5_base_covid --export=ALL,VIRUS=covid,METHODS="$METH" scripts/slurm_table5_baselines.sh
submit --job-name=t5_base_h1 --export=ALL,VIRUS=h1n1,METHODS="$METH" scripts/slurm_table5_baselines.sh

# Table 8 COVID ablation evals (gen flags on v5)
for A in no_bridge no_tree_ctx no_seq_branch no_bl_eval full; do
  submit --job-name=t8_$A --export=ALL,ABLATE=$A scripts/slurm_table8_ablations.sh
done

# Table 8 retrain: no BL head
submit --job-name=t8_nobl --time=24:00:00 --export=ALL,ABLATE=train_no_bl scripts/slurm_table8_ablations.sh

# Optional COVID baselines table (Tree-KL etc.) if time
# submit --job-name=t8_base_tab --export=ALL,DATA=data/covid/test,TRAIN_DATA=data/covid/train,OUT=benchmarks/results/results_baselines_covid.csv \
#   scripts/slurm_baselines.sh checkpoints/covid_v5_mutrec/best.pt

echo "==== job ids ===="
printf '%s\n' "${JOBS[@]}"
mkdir -p benchmarks/results/tables
printf '%s\n' "${JOBS[@]}" > benchmarks/results/tables/table5_table8_job_ids.txt
python3 - <<'PY'
import json
from pathlib import Path
ids = Path("benchmarks/results/tables/table5_table8_job_ids.txt").read_text().split()
Path("benchmarks/results/tables/table8_job_ids.json").write_text(json.dumps({
  "submitted": ids,
  "note": "See squeue / sacct for mapping by job name t5_* / t8_*",
}, indent=2))
print("wrote table8_job_ids.json", ids)
PY
squeue -u nnori | head -40
REMOTE
