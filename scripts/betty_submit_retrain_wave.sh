#!/bin/bash
# Master submit for Retrain + VaxSeer + Ab SHM wave (Betty login03).
# Sync code first, then:
#   bash scripts/betty_submit_retrain_wave.sh
#
# Order: (1) viral eval sweep  (2) viral retrain  (3) Ab v2 SHM train+rod
#         (4) H3 tree viz + variant path  (5) VaxSeer setup (non-sbatch)
# Does NOT overwrite paper ckpts.

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs benchmarks/results/tables

echo "=== 1. Viral eval-only sweep ==="
bash scripts/betty_submit_viral_eval_sweep.sh

echo "=== 2. Viral mutlin retrain (β=0 and 0.25) ==="
bash scripts/betty_submit_viral_retrain.sh

echo "=== 3. Ab OAS v2 SHM train + afterok Rod.82 ==="
TRAIN_J=$(sbatch --parsable --qos=mig-max scripts/slurm_ab_oas_v2_shm_train.sh)
echo "ab_oas_v2_shm train -> $TRAIN_J"
ROD_J=$(sbatch --parsable --dependency=afterok:${TRAIN_J} scripts/slurm_ab_oas_v2_shm_rod82.sh)
echo "ab_oas_v2_shm rod82 -> $ROD_J"
python - <<PY
import json
from pathlib import Path
p = Path("benchmarks/results/tables/ab_oas_v2_shm_job_ids.json")
p.write_text(json.dumps({"train": int("$TRAIN_J"), "rod82": int("$ROD_J")}, indent=2) + "\n")
print("Wrote", p)
PY

echo "=== 4. H3 tree viz (for variant-path figure) ==="
VIZ_J=$(sbatch --parsable --qos=mig-max scripts/slurm_h3n2_tree_viz.sh)
echo "h3n2_tree_viz -> $VIZ_J"
sbatch --parsable --dependency=afterok:${VIZ_J} --partition=genoa-std-mem \
  --cpus-per-task=4 --mem-per-cpu=5632M --time=2:00:00 \
  --job-name=h3_varpath --output=logs/h3_variant_path_%j.log \
  --wrap='cd ~/DiscreteTreeFlows; /vast/home/n/nnori/.conda/envs/treesbm/bin/python scripts/build_flu_variant_path_figure.py --virus h3n2 --clade-name 3C.2a1b.2a.2'

echo "=== 5. Ab CDR/FWR entropy audit (CPU) ==="
sbatch --parsable --partition=genoa-std-mem --cpus-per-task=4 --mem-per-cpu=5632M --time=4:00:00 \
  --job-name=ab_entropy --output=logs/ab_entropy_audit_%j.log \
  --wrap='cd ~/DiscreteTreeFlows; export PYTHONPATH=$HOME/DiscreteTreeFlows:$PYTHONPATH; /vast/home/n/nnori/.conda/envs/treesbm/bin/python scripts/ab_cdr_fwr_entropy_audit.py'

echo "=== 6. VaxSeer: clone done under \$LABHOME/third_party/vaxseer; lm download needs Dropbox token ==="
echo "  cd \$LABHOME/third_party/vaxseer"
echo "  python download_models_from_dropbox.py --task lm --year 2018 --subtype a_h3n2 --output_dir runs"
echo "  See benchmarks/results/tables/vaxseer_dominance.md"

echo "=== wave submit done ==="
squeue -u "$USER" | head -40
