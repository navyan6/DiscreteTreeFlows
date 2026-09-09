#!/bin/bash
#SBATCH --job-name=eval_h1n1_lit
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/eval_h1n1_lit_%j.log
#SBATCH --error=logs/eval_h1n1_lit_%j.log
#
# Full-metric H1N1 lit-hotspot enrichment (H1 nmicrobiol mask + EVEscape flu_h1).
# HA head/RBS region bands are H3-numbered — disabled here (pathogen separation).
#   EVAL_MRS="0.3 0.5" OUT_SUFFIX=_fullmetrics \
#     sbatch --qos=mig-max scripts/_run_eval_h1n1_lit.sh
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
export PYTHONUNBUFFERED=1
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
cd ~/DiscreteTreeFlows
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
FLU_MASK="${FLU_HOTSPOT_MASK:-results/flu_mutfreq_vs_lit/mut_hotspot_mask_h1_nmicrobiol_lit.pt}"
CKPT_NAME="${CKPT_NAME:-h1n1_v2_lit_hotspot}"
CKPT="checkpoints/${CKPT_NAME}/best.pt"
EVAL_MRS="${EVAL_MRS:-0.3 0.5}"
OUT_SUFFIX="${OUT_SUFFIX:-_fullmetrics}"
SKIP_EXISTING="${SKIP_EXISTING:-0}"
NSTEPS="${NSTEPS:-100}"
MAX_TREES="${MAX_TREES:-20}"
EVESCAPE="${EVESCAPE:-data/evescape_h1n1_ha.pt}"
DATA="${DATA:-data/h1n1/test}"

test -f "$CKPT" || { echo "ERROR: missing $CKPT"; exit 1; }
test -f "$FLU_MASK" || { echo "ERROR: missing $FLU_MASK"; exit 1; }
test -f "$EVESCAPE" || { echo "ERROR: missing $EVESCAPE — prepare via prepare_evescape.py from flu_h1_evescape.csv"; exit 1; }
echo "Start: $(date)  ckpt=$CKPT  data=$DATA  mrs=[$EVAL_MRS]  OUT_SUFFIX=$OUT_SUFFIX"
ls -lah "$CKPT" "$FLU_MASK" "$EVESCAPE"
echo "Metrics: mut_recovery, site_recall, aa_acc_given_hit, site_precision,"
echo "  cons_retention, identity, dist_to_root, gt_dist_to_root,"
echo "  mut_recovery_any_descendant, site_recall_any_descendant, mut_site_recall_path_union,"
echo "  flu_hotspot_mut_frac / lit_hotspot_mut_frac,"
echo "  evescape_mean_antigenic_muts / gt_evescape_mean_antigenic_muts / random_baseline_evescape_antigenic,"
echo "  evescape_mean_all_scored_muts (= model_evescape) companion,"
echo "  (mean of per-mut EVEscape at antigenic sites — not a strain product)"
echo "HA head/RBS region bands: OFF (H3-numbered; H1 uses lit mask only)"

for MRS in $EVAL_MRS; do
  OUT="checkpoints/eval_enrichment_${CKPT_NAME}_mrs${MRS}${OUT_SUFFIX}.json"
  if [[ "$SKIP_EXISTING" == "1" && -f "$OUT" ]]; then
    echo "SKIP enrichment mrs=$MRS (exists: $OUT)"
    continue
  fi
  echo "==== mrs=$MRS -> $OUT ===="
  "$PYTHON" -u scripts/eval_evescape_enrichment.py \
    --checkpoint "$CKPT" \
    --data "$DATA" \
    --max-seq-len 566 \
    --lit-hotspot-mask "$FLU_MASK" \
    --no-ha-region-metrics \
    --evescape "$EVESCAPE" \
    --mutation-rate-scale "$MRS" \
    --n-steps "$NSTEPS" \
    --max-trees "$MAX_TREES" \
    --out "$OUT"
done
echo "Done: $(date)"
