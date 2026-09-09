#!/bin/bash
#SBATCH --job-name=eval_h3n2_lit
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/eval_h3n2_lit_%j.log
#SBATCH --error=logs/eval_h3n2_lit_%j.log
#
# Full-metric H3N2 lit-hotspot enrichment (nmicrobiol mask + HA head/RBS
# region metrics + dist_to_root + optional EVEscape if tensor present).
#   EVAL_MRS="0.3 0.5" OUT_SUFFIX=_fullmetrics SKIP_EXISTING=0 \
#     sbatch --qos=mig-max scripts/_run_eval_h3n2_lit.sh
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
export PYTHONUNBUFFERED=1
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
cd ~/DiscreteTreeFlows
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
FLU_MASK="${FLU_HOTSPOT_MASK:-results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt}"
CKPT_NAME="${CKPT_NAME:-h3n2_v3_lit_hotspot}"
CKPT="checkpoints/${CKPT_NAME}/best.pt"
EVAL_MRS="${EVAL_MRS:-0.3 0.5}"
OUT_SUFFIX="${OUT_SUFFIX:-_fullmetrics}"
SKIP_EXISTING="${SKIP_EXISTING:-0}"
NSTEPS="${NSTEPS:-100}"
MAX_TREES="${MAX_TREES:-20}"
# Official EVEscape release has flu_h1 only — H3 tensor is optional.
EVESCAPE="${EVESCAPE:-data/evescape_h3n2_ha.pt}"

test -f "$CKPT" || { echo "ERROR: missing $CKPT"; exit 1; }
test -f "$FLU_MASK" || { echo "ERROR: missing $FLU_MASK"; exit 1; }
echo "Start: $(date)  ckpt=$CKPT  mrs=[$EVAL_MRS]  OUT_SUFFIX=$OUT_SUFFIX  SKIP_EXISTING=$SKIP_EXISTING"
ls -lah "$CKPT" "$FLU_MASK"
echo "Metrics (always): mut_recovery, site_recall, aa_acc_given_hit, site_precision,"
echo "  cons_retention, identity, dist_to_root, gt_dist_to_root,"
echo "  mut_recovery_any_descendant, site_recall_any_descendant, mut_site_recall_path_union,"
echo "  flu_hotspot_mut_frac / lit_hotspot_mut_frac,"
echo "  ha_head_{mut_frac,site_recall,any_mut}, ha_rbs_{mut_frac,site_recall,any_mut}"
EV_ARGS=()
if [[ -f "$EVESCAPE" ]]; then
  EV_ARGS+=(--evescape "$EVESCAPE")
  echo "EVEscape: ON ($EVESCAPE)"
else
  echo "EVEscape: OFF (missing $EVESCAPE — official release has flu_h1 only; H3 optional)"
fi

for MRS in $EVAL_MRS; do
  OUT="checkpoints/eval_enrichment_${CKPT_NAME}_mrs${MRS}${OUT_SUFFIX}.json"
  if [[ "$SKIP_EXISTING" == "1" && -f "$OUT" ]]; then
    echo "SKIP enrichment mrs=$MRS (exists: $OUT)"
    continue
  fi
  echo "==== mrs=$MRS -> $OUT ===="
  "$PYTHON" -u scripts/eval_evescape_enrichment.py \
    --checkpoint "$CKPT" \
    --data data/h3n2/test \
    --max-seq-len 566 \
    --lit-hotspot-mask "$FLU_MASK" \
    --ha-region-metrics \
    --mutation-rate-scale "$MRS" \
    --n-steps "$NSTEPS" \
    --max-trees "$MAX_TREES" \
    --out "$OUT" \
    "${EV_ARGS[@]}"
done
echo "Done: $(date)"
echo "Outputs: checkpoints/eval_enrichment_${CKPT_NAME}_mrs*{${OUT_SUFFIX}}.json"
