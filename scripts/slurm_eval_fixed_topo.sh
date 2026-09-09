#!/bin/bash
#SBATCH --job-name=ftopo_seq
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/eval_fixed_topo_%j.log
#SBATCH --error=logs/eval_fixed_topo_%j.log
#
# Fixed-topology dest bake-off. Read-only on --checkpoint.
#   VIRUS=covid CKPT=checkpoints/covid_v5_mutrec/best.pt \
#     sbatch --qos=mig-max scripts/slurm_eval_fixed_topo.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints/viral_eval_sweep

VIRUS="${VIRUS:-covid}"
MRS="${MRS:-1.0}"
NT_FIRE="${NT_FIRE:-1.0}"
MAX_TREES="${MAX_TREES:-20}"
METHODS="${METHODS:-neutral,codon_gy94,treesbm,treesbm_ntfire,codon_gy94_ntfire}"

case "$VIRUS" in
  covid)
    DATA="${DATA:-data/covid/test}"
    MAX_SEQ_LEN="${MAX_SEQ_LEN:-1280}"
    CKPT="${CKPT:-checkpoints/covid_v5_mutrec/best.pt}"
    ;;
  h1n1)
    DATA="${DATA:-data/h1n1/test}"
    MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"
    CKPT="${CKPT:-checkpoints/h1n1_v2_lit_hotspot/best.pt}"
    ;;
  h3n2)
    DATA="${DATA:-data/h3n2/test}"
    MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"
    CKPT="${CKPT:-checkpoints/h3n2_v3_lit_hotspot/best.pt}"
    ;;
  hiv_geo)
    DATA="${DATA:-data/hiv_geo/test}"
    MAX_SEQ_LEN="${MAX_SEQ_LEN:-900}"
    CKPT="${CKPT:-checkpoints/hiv_geo_v1/best.pt}"
    ;;
  bdbv)
    DATA="${DATA:-data/bdbv_temporal/test}"
    if [[ -f results/bdbv_l_conservation/window_config.json ]]; then
      MAX_SEQ_LEN="${MAX_SEQ_LEN:-$($PYTHON -c "import json; print(json.load(open('results/bdbv_l_conservation/window_config.json'))['max_seq_len'])")}"
    else
      MAX_SEQ_LEN="${MAX_SEQ_LEN:-900}"
    fi
    CKPT="${CKPT:-checkpoints/bdbv_v1_mutrec/best.pt}"
    ;;
  *)
    echo "Unknown VIRUS=$VIRUS" >&2; exit 2 ;;
esac

OUT="${OUT:-checkpoints/viral_eval_sweep/${VIRUS}_fixed_topo_mrs${MRS}.json}"
echo "virus=$VIRUS ckpt=$CKPT data=$DATA out=$OUT"
test -f "$CKPT" || { echo "missing ckpt $CKPT" >&2; exit 1; }
test -d "$DATA" || { echo "missing data $DATA" >&2; exit 1; }

$PYTHON -u scripts/eval_fixed_topo_seq.py \
  --data "$DATA" \
  --max-seq-len "$MAX_SEQ_LEN" \
  --checkpoint "$CKPT" \
  --methods "$METHODS" \
  --mutation-rate-scale "$MRS" \
  --nt-fire-scale "$NT_FIRE" \
  --max-trees "$MAX_TREES" \
  --out "$OUT"

echo "Done $(date -Is) $OUT"
