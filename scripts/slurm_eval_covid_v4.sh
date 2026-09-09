#!/bin/bash
#SBATCH --job-name=eval_covid_v4
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/eval_covid_v4_%j.log
#SBATCH --error=logs/eval_covid_v4_%j.log
#
# Eval covid_v4_mutrec (pos-emb + λ_mut=8 + λ_semi) at mrs ∈ {0.3, 0.5, 1.0}.
# Safe to run while training job still writes later checkpoints — reads best.pt.
# SKIP_EXISTING=1 avoids duplicating Agent G enrichment JSONs.
#
#   sbatch --qos=mig-max scripts/slurm_eval_covid_v4.sh
#   EVAL_MRS="1.0" SKIP_EXISTING=1 sbatch --qos=mig-max scripts/slurm_eval_covid_v4.sh

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export PYTHONUNBUFFERED=1
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

CKPT_NAME="${CKPT_NAME:-covid_v4_mutrec}"
CKPT="checkpoints/${CKPT_NAME}/best.pt"
if [[ ! -f "$CKPT" && -f "$LABHOME/checkpoints_backup/${CKPT_NAME}/best.pt" ]]; then
    mkdir -p "checkpoints/${CKPT_NAME}"
    cp -n "$LABHOME/checkpoints_backup/${CKPT_NAME}/best.pt" "$CKPT"
fi
if [[ ! -f "$CKPT" ]]; then
    echo "ERROR: $CKPT missing — wait for covid_mutrec training or copy from LABHOME backup"
    exit 1
fi

EVESCAPE=data/covid/evescape_spike_rbd.pt
EVAL_MRS="${EVAL_MRS:-0.3 0.5 1.0}"
NSTEPS="${NSTEPS:-100}"
MAX_TREES="${MAX_TREES:-20}"
SKIP_EXISTING="${SKIP_EXISTING:-1}"
OUT_SUFFIX="${OUT_SUFFIX:-}"  # e.g. _regions → ..._mrs0.3_regions.json

echo "Start: $(date)  ckpt=$CKPT  mrs=[$EVAL_MRS]  SKIP_EXISTING=$SKIP_EXISTING  OUT_SUFFIX=$OUT_SUFFIX"
ls -lah "$CKPT"

for MRS in $EVAL_MRS; do
    OUT="checkpoints/eval_enrichment_${CKPT_NAME}_mrs${MRS}${OUT_SUFFIX}.json"
    if [[ "$SKIP_EXISTING" == "1" && -f "$OUT" ]]; then
        echo "SKIP enrichment mrs=$MRS (exists: $OUT)"
        continue
    fi
    echo ""
    echo "############################################################"
    echo "# enrichment  $CKPT_NAME  mrs=$MRS"
    echo "############################################################"
    PMC_MASK="${PMC_HOTSPOT_MASK:-results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt}"
    PMC_ARGS=()
    if [[ -n "$PMC_MASK" && -f "$PMC_MASK" ]]; then
        # Explicit COVID-only path (never flu masks).
        PMC_ARGS+=(--lit-hotspot-mask "$PMC_MASK")
    else
        PMC_ARGS+=(--lit-hotspot-mask "")
    fi
    $PYTHON -u scripts/eval_evescape_enrichment.py \
        --checkpoint "$CKPT" \
        --data data/covid/test \
        --max-seq-len 1280 \
        --evescape "$EVESCAPE" \
        --mutation-rate-scale "$MRS" \
        --n-steps "$NSTEPS" \
        --max-trees "$MAX_TREES" \
        --out "$OUT" \
        "${PMC_ARGS[@]}"
    if [[ "$MRS" == "0.3" && -z "$OUT_SUFFIX" && -f "$OUT" ]]; then
        cp -f "$OUT" "checkpoints/eval_enrichment_${CKPT_NAME}.json"
    fi
done

echo "Done: $(date)"
