#!/bin/bash
#SBATCH --job-name=eval_covid
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/eval_covid_%j.log
#SBATCH --error=logs/eval_covid_%j.log
#
# Covid EVEscape enrichment + mut/cons recovery (Agent G Wave-2).
# Default QOS mig-max avoids MaxGRESPerAccount on qos=mig.
#
# Priority: covid_v3_cons at mrs ∈ {0.3, 0.5, 1.0} (Wave-1 only finished mrs=0.1).
# Optional: EVAL_EXTRA_CKPTS="covid_v1 covid_v2_entropy" to also score older ckpts.
# Optional: EVAL_MRS="0.3 0.5 1.0"  EVAL_MAX_TREES=20  SITE_SOFTMAX=1
# Skips existing enrichment JSONs (SKIP_EXISTING=1) so Agent G jobs are not duplicated.
#
# EVEscape tensor required at: data/covid/evescape_spike_rbd.pt

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
export TRANSFORMERS_CACHE=${TRANSFORMERS_CACHE:-$HF_HOME}
export TORCH_HOME=${TORCH_HOME:-$LABHOME/torch_cache}

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

EVESCAPE=data/covid/evescape_spike_rbd.pt
if [[ ! -f "$EVESCAPE" ]]; then
    echo "ERROR: EVEscape tensor missing at $EVESCAPE"
    exit 1
fi
echo "EVEscape OK: $EVESCAPE ($(du -h "$EVESCAPE" | awk '{print $1}'))"

resolve_ckpt() {
    local name="$1"
    if [[ -f "checkpoints/${name}/best.pt" ]]; then
        echo "checkpoints/${name}/best.pt"
    elif [[ -f "$LABHOME/checkpoints_backup/${name}/best.pt" ]]; then
        mkdir -p "checkpoints/${name}"
        cp -n "$LABHOME/checkpoints_backup/${name}/best.pt" "checkpoints/${name}/best.pt"
        echo "checkpoints/${name}/best.pt"
    else
        echo ""
    fi
}

EVAL_MRS="${EVAL_MRS:-0.3 0.5 1.0}"
EVAL_MAX_TREES="${EVAL_MAX_TREES:-20}"
NSTEPS="${NSTEPS:-100}"
SKIP_EXISTING="${SKIP_EXISTING:-1}"
SITE_SOFTMAX="${SITE_SOFTMAX:-0}"
SITE_TEMPERATURE="${SITE_TEMPERATURE:-1.0}"
# Keep Wave-1 timeout mode away: do not load v1/v2 unless explicitly requested.
CKPTS="${EVAL_CKPTS:-covid_v3_cons}"
EXTRA="${EVAL_EXTRA_CKPTS:-}"

SOFTMAX_ARGS=()
SOFTMAX_TAG=""
if [[ "$SITE_SOFTMAX" == "1" || "$SITE_SOFTMAX" == "true" ]]; then
    SOFTMAX_ARGS+=(--site-softmax-sample --site-temperature "$SITE_TEMPERATURE")
    SOFTMAX_TAG="_sitesm"
fi

echo "Start: $(date)  LABHOME=$LABHOME  mrs=[$EVAL_MRS] ckpts=[$CKPTS $EXTRA]"
echo "SKIP_EXISTING=$SKIP_EXISTING  SITE_SOFTMAX=$SITE_SOFTMAX  nsteps=$NSTEPS"

for ckpt in $CKPTS $EXTRA; do
    CKPT_PATH=$(resolve_ckpt "$ckpt")
    if [[ -z "$CKPT_PATH" ]]; then
        echo "SKIP $ckpt — checkpoint not found locally or in LABHOME backup"
        continue
    fi
    for MRS in $EVAL_MRS; do
        OUT="checkpoints/eval_enrichment_${ckpt}_mrs${MRS}${SOFTMAX_TAG}.json"
        # Also write legacy name for mrs=0.3 (playbook / table consumers).
        if [[ "$MRS" == "0.3" && -z "$SOFTMAX_TAG" ]]; then
            OUT_LEGACY="checkpoints/eval_enrichment_${ckpt}.json"
        else
            OUT_LEGACY=""
        fi
        if [[ "$SKIP_EXISTING" == "1" && -f "$OUT" ]]; then
            echo "SKIP enrichment $ckpt mrs=$MRS (exists: $OUT)"
            continue
        fi
        echo ""
        echo "############################################################"
        echo "# enrichment  $ckpt  mrs=$MRS  ($CKPT_PATH)  sitesm=$SITE_SOFTMAX"
        echo "############################################################"
        PMC_MASK="${PMC_HOTSPOT_MASK:-results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt}"
        PMC_ARGS=()
        if [[ -n "$PMC_MASK" && -f "$PMC_MASK" ]]; then
            # Explicit COVID-only path (never flu masks). Prefer --lit-hotspot-mask.
            PMC_ARGS+=(--lit-hotspot-mask "$PMC_MASK")
        else
            PMC_ARGS+=(--lit-hotspot-mask "")
        fi
        $PYTHON scripts/eval_evescape_enrichment.py \
            --checkpoint "$CKPT_PATH" \
            --data data/covid/test \
            --max-seq-len 1280 \
            --evescape "$EVESCAPE" \
            --mutation-rate-scale "$MRS" \
            --n-steps "$NSTEPS" \
            --max-trees "$EVAL_MAX_TREES" \
            --out "$OUT" \
            "${PMC_ARGS[@]}" \
            "${SOFTMAX_ARGS[@]}"
        if [[ -n "$OUT_LEGACY" && -f "$OUT" ]]; then
            cp -f "$OUT" "$OUT_LEGACY"
            echo "Also wrote $OUT_LEGACY"
        fi
    done
done

# Smoke deep metrics on covid_v3_cons @ mrs=0.3 for 1–2 test groups.
CKPT_V3=$(resolve_ckpt covid_v3_cons)
if [[ -n "$CKPT_V3" ]]; then
    mapfile -t GROUPS < <(ls data/covid/test/group_*_rooted.nwk 2>/dev/null \
        | sed -E 's/.*group_0*([0-9]+)_rooted\.nwk/\1/' | head -2)
    if [[ ${#GROUPS[@]} -eq 0 ]]; then
        echo "SKIP smoke eval_single_tree — no group_*_rooted.nwk under data/covid/test"
    else
        for g in "${GROUPS[@]}"; do
            echo ""
            echo "############################################################"
            echo "# smoke eval_single_tree  covid_v3_cons  group=$g  mrs=0.3"
            echo "############################################################"
            $PYTHON scripts/eval_single_tree.py \
                --checkpoint "$CKPT_V3" \
                --data data/covid/test \
                --group "$g" \
                --max-seq-len 1280 \
                --max-leaves 300 \
                --mutation-rate-scale 0.3 \
                --n-steps 100
        done
    fi
fi

echo ""
echo "Done: $(date)"
echo "JSON outputs: checkpoints/eval_enrichment_<ckpt>_mrs<mrs>.json"
