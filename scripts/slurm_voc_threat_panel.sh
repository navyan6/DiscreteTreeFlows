#!/bin/bash
#SBATCH --job-name=voc_threat
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/voc_threat_%j.log
#SBATCH --error=logs/voc_threat_%j.log
#
# Generate TreeSBM trees for one VOC threat case (via eval_single_tree), then
# score VOC-bundle metrics with eval_voc_threat_recovery.py.
#
# Required env:
#   CASE_DIR=results/voc_threat_panel/cases/Gamma_g003_test
#   VOC_ID=Gamma
#   DATA=data/covid/test
#   GROUP=3
#
# Optional:
#   CKPT_NAME=covid_v7_pmc_hotspot  MRS=0.5  MAX_LEAVES=250  N_STEPS=160
#   SEED=42  BRANCH_SCALE=6.0  REF_TIP=MZ397166.1  TOPK=10
#
# Usage:
#   CASE_DIR=... VOC_ID=Gamma DATA=data/covid/test GROUP=3 \
#     sbatch scripts/slurm_voc_threat_panel.sh

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
export TRANSFORMERS_CACHE=${TRANSFORMERS_CACHE:-$HF_HOME}
export TORCH_HOME=${TORCH_HOME:-$LABHOME/torch_cache}

cd ~/DiscreteTreeFlows
mkdir -p logs results/voc_threat_panel checkpoints

CASE_DIR="${CASE_DIR:?set CASE_DIR}"
VOC_ID="${VOC_ID:?set VOC_ID}"
DATA="${DATA:?set DATA (e.g. data/covid/test)}"
GROUP="${GROUP:?set GROUP (integer)}"
CKPT_NAME="${CKPT_NAME:-covid_v7_pmc_hotspot}"
CKPT="checkpoints/${CKPT_NAME}/best.pt"
FALLBACK="checkpoints/covid_v5_mutrec/best.pt"
MRS="${MRS:-0.5}"
MAX_LEAVES="${MAX_LEAVES:-250}"
N_STEPS="${N_STEPS:-160}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-1280}"
BRANCH_SCALE="${BRANCH_SCALE:-6.0}"
SEED="${SEED:-42}"
REF_TIP="${REF_TIP:-}"
TOPK="${TOPK:-10}"

if [[ ! -f "$CKPT" ]]; then
  echo "WARN: $CKPT missing — falling back to $FALLBACK"
  CKPT="$FALLBACK"
  CKPT_NAME="covid_v5_mutrec"
fi
if [[ ! -f "$CKPT" ]]; then
  echo "ERROR: no usable checkpoint"
  exit 1
fi

mkdir -p "$CASE_DIR"
OBS_FA="$DATA/group_$(printf '%03d' "$GROUP")_anc_aa.fasta"
OBS_NWK="$DATA/group_$(printf '%03d' "$GROUP")_rooted.nwk"
if [[ ! -f "$OBS_FA" ]]; then
  echo "ERROR: missing observed fasta $OBS_FA"
  exit 1
fi
cp -f "$OBS_FA" "$CASE_DIR/observed_anc_aa.fasta"
[[ -f "$OBS_NWK" ]] && cp -f "$OBS_NWK" "$CASE_DIR/observed.nwk"

echo "Start: $(date)  VOC=$VOC_ID  DATA=$DATA  GROUP=$GROUP  CKPT=$CKPT_NAME"
echo "MRS=$MRS MAX_LEAVES=$MAX_LEAVES N_STEPS=$N_STEPS SEED=$SEED"

$PYTHON scripts/eval_single_tree.py \
    --checkpoint "$CKPT" \
    --data "$DATA" \
    --group "$GROUP" \
    --max-seq-len "$MAX_SEQ_LEN" \
    --max-leaves "$MAX_LEAVES" \
    --mutation-rate-scale "$MRS" \
    --n-steps "$N_STEPS" \
    --branch-rate-scale "$BRANCH_SCALE" \
    --seed "$SEED"

SRC_NWK="checkpoints/gen_group${GROUP}.nwk"
SRC_FA="checkpoints/gen_group${GROUP}.fasta"
if [[ ! -f "$SRC_NWK" || ! -f "$SRC_FA" ]]; then
  echo "ERROR: expected $SRC_NWK and $SRC_FA after eval_single_tree"
  ls -la checkpoints/gen_group${GROUP}.* 2>/dev/null || true
  exit 1
fi

GEN_NWK="$CASE_DIR/generated.nwk"
GEN_FA="$CASE_DIR/generated.fasta"
cp -f "$SRC_NWK" "$GEN_NWK"
cp -f "$SRC_FA" "$GEN_FA"
cp -f "$SRC_NWK" "$CASE_DIR/generated_s${SEED}.nwk"
cp -f "$SRC_FA" "$CASE_DIR/generated_s${SEED}.fasta"

EVAL_OUT="$CASE_DIR/voc_eval_${VOC_ID}.json"
REF_ARGS=()
if [[ -n "$REF_TIP" ]]; then
  REF_ARGS+=(--ref-tip "$REF_TIP")
fi

$PYTHON scripts/eval_voc_threat_recovery.py \
  --voc "$VOC_ID" \
  --gen-fasta "$GEN_FA" \
  --obs-fasta "$CASE_DIR/observed_anc_aa.fasta" \
  --gen-nwk "$GEN_NWK" \
  --obs-nwk "$CASE_DIR/observed.nwk" \
  --topk "$TOPK" \
  --evescape data/covid/evescape_spike_rbd.pt \
  --out "$EVAL_OUT" \
  "${REF_ARGS[@]}"

META="$CASE_DIR/run_meta.json"
$PYTHON - <<PY
import json
from pathlib import Path
Path("$META").write_text(json.dumps({
    "voc_id": "$VOC_ID",
    "data": "$DATA",
    "group": int("$GROUP"),
    "ckpt": "$CKPT",
    "ckpt_name": "$CKPT_NAME",
    "mrs": float("$MRS"),
    "max_leaves": int("$MAX_LEAVES"),
    "n_steps": int("$N_STEPS"),
    "branch_scale": float("$BRANCH_SCALE"),
    "seed": int("$SEED"),
    "slurm_job": "${SLURM_JOB_ID:-local}",
    "eval_out": "$EVAL_OUT",
}, indent=2) + "\n")
PY

echo "Done: $(date)  -> $EVAL_OUT"
$PYTHON - <<PY
import json
d=json.load(open("$EVAL_OUT"))
g=d["gen_vs_score_root"]
print("exact", round(g["voc_exact_mut_recall"],3),
      "acq", round(g["voc_acquired_mut_recall"],3),
      "bundle_any", g["voc_bundle_any"],
      "topk_union", d["topk"]["union_sig"])
PY
