#!/bin/bash
#SBATCH --job-name=c1_ts_nll
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/c1_treesbm_plmnll_%j.log
#SBATCH --error=logs/c1_treesbm_plmnll_%j.log
#
# C.1 TreeSBM free-gen enrich with pLM NLL (default-on in eval_evescape_enrichment).
#
#   CKPT=checkpoints/h3n2_v3_lit_hotspot/best.pt DATA=data/h3n2/test MAX_SEQ_LEN=566 \
#     OUT=checkpoints/eval_enrichment_h3n2_v3_lit_hotspot_mrs0.5_plmnll.json \
#     LIT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt \
#     sbatch --qos=mig-max scripts/slurm_c1_treesbm_plmnll.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

CKPT="${CKPT:-checkpoints/h3n2_v3_lit_hotspot/best.pt}"
DATA="${DATA:-data/h3n2/test}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"
MRS="${MRS:-0.5}"
NSTEPS="${NSTEPS:-100}"
MAX_TREES="${MAX_TREES:-20}"
OUT="${OUT:-checkpoints/eval_enrichment_h3n2_v3_lit_hotspot_mrs0.5_plmnll.json}"
LIT_MASK="${LIT_MASK:-}"

resolve_ckpt() {
  local c="$1"
  if [[ -f "$c" ]]; then echo "$c"; return; fi
  local name; name=$(basename "$(dirname "$c")")
  local backup="$LABHOME/checkpoints_backup/${name}/best.pt"
  if [[ -f "$backup" ]]; then
    mkdir -p "$(dirname "$c")"
    cp -n "$backup" "$c" 2>/dev/null || cp "$backup" "$c"
    echo "$c"; return
  fi
  echo "$c"
}

CKPT_PATH=$(resolve_ckpt "$CKPT")
MASK_ARGS=()
if [[ -n "$LIT_MASK" && -f "$LIT_MASK" ]]; then
  MASK_ARGS=(--lit-hotspot-mask "$LIT_MASK")
fi

echo "Start: $(date) ckpt=$CKPT_PATH data=$DATA L=$MAX_SEQ_LEN mrs=$MRS out=$OUT"
$PYTHON scripts/eval_evescape_enrichment.py \
  --checkpoint "$CKPT_PATH" \
  --data "$DATA" \
  --max-seq-len "$MAX_SEQ_LEN" \
  --mutation-rate-scale "$MRS" \
  --n-steps "$NSTEPS" \
  --max-trees "$MAX_TREES" \
  --out "$OUT" \
  "${MASK_ARGS[@]}"

echo "Done: $(date) -> $OUT"
python3 - <<PY
import json
from pathlib import Path
p = Path("$OUT")
d = json.loads(p.read_text())
s = d.get("summary", d)
print("plm_nll=", s.get("plm_nll"))
PY
