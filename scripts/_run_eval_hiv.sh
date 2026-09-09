#!/bin/bash
#SBATCH --job-name=eval_hiv
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/eval_hiv_%j.log
#SBATCH --error=logs/eval_hiv_%j.log
#
# HIV-1 Env full-metric enrichment (V1–V5 lit hotspot + optional EVEscape).
# Mirrors scripts/_run_eval_h1n1_lit.sh / _run_eval_h3n2_lit.sh.
#
# Usage:
#   CKPT_NAME=hiv_temporal_v1 DATA=data/hiv_temporal/test \
#     sbatch --qos=mig-max scripts/_run_eval_hiv.sh
#   CKPT_NAME=hiv_geo_v1 DATA=data/hiv_geo/test \
#     sbatch --qos=mig-max --dependency=afterok:<train_job> scripts/_run_eval_hiv.sh
#
# Locked:
#   max_seq_len=900 (Env ~856–898; do NOT use flu 566)
#   HA head/RBS and Spike domain metrics OFF
#   hotspot = results/hiv_env_mask/mut_hotspot_mask_hiv_env_v1v5.pt
#   EVEscape ON only if data/evescape_hiv_env.pt exists (else prepare from
#   EVEscape/results/summaries_with_scores/hiv_env_evescape.csv; else OFF)
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
export PYTHONUNBUFFERED=1
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
cd ~/DiscreteTreeFlows
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

CKPT_NAME="${CKPT_NAME:?set CKPT_NAME=hiv_temporal_v1 or hiv_geo_v1}"
CKPT="checkpoints/${CKPT_NAME}/best.pt"
DATA="${DATA:?set DATA=data/hiv_temporal/test or data/hiv_geo/test}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-900}"
HIV_MASK="${HIV_HOTSPOT_MASK:-results/hiv_env_mask/mut_hotspot_mask_hiv_env_v1v5.pt}"
EVAL_MRS="${EVAL_MRS:-0.3 0.5}"
OUT_SUFFIX="${OUT_SUFFIX:-_fullmetrics}"
SKIP_EXISTING="${SKIP_EXISTING:-0}"
NSTEPS="${NSTEPS:-100}"
MAX_TREES="${MAX_TREES:-20}"
EVESCAPE="${EVESCAPE:-data/evescape_hiv_env.pt}"
EVESCAPE_CSV="${EVESCAPE_CSV:-EVEscape/results/summaries_with_scores/hiv_env_evescape.csv}"

test -f "$CKPT" || { echo "ERROR: missing $CKPT"; exit 1; }
test -d "$DATA" || { echo "ERROR: missing $DATA"; exit 1; }
test -f "$HIV_MASK" || { echo "ERROR: missing $HIV_MASK"; exit 1; }

# Refuse flu/COVID masks and flu-length eval
case "$HIV_MASK" in
  *covid_mutfreq*|*pmc_lit*|*flu_mutfreq*|*nmicrobiol*|*ab_cdr*)
    echo "REFUSING non-HIV hotspot mask: $HIV_MASK" >&2
    exit 2
    ;;
esac
if [[ "$MAX_SEQ_LEN" -le 600 ]]; then
  echo "REFUSING max_seq_len=$MAX_SEQ_LEN (flu HA). HIV Env must be ~900." >&2
  exit 2
fi

echo "Start: $(date)  ckpt=$CKPT  data=$DATA  mrs=[$EVAL_MRS]  L=$MAX_SEQ_LEN"
ls -lah "$CKPT" "$HIV_MASK"
echo "Metrics: mut_recovery, site_recall, aa_acc_given_hit, site_precision,"
echo "  cons_retention, identity, dist_to_root, gt_dist_to_root,"
echo "  mut_recovery_any_descendant, site_recall_any_descendant, mut_site_recall_path_union,"
echo "  lit_hotspot_mut_frac (HXB2 Env V1–V5),"
echo "  evescape_mean_antigenic_muts if HIV Env tensor present"
echo "HA head/RBS + Spike domain bands: OFF (HIV Env, not flu/COVID)"

# Prepare HIV EVEscape tensor if CSV exists and .pt is missing.
if [[ ! -f "$EVESCAPE" && -f "$EVESCAPE_CSV" ]]; then
  echo "Preparing HIV EVEscape tensor from $EVESCAPE_CSV"
  REF_FA="data/_hiv_evescape_ref.fa"
  "$PYTHON" - <<PY
from pathlib import Path
from Bio import SeqIO
best_id, best = "hiv_env_ref", ""
for f in Path("$DATA").glob("group_*_anc_aa.fasta"):
    for r in SeqIO.parse(str(f), "fasta"):
        s = str(r.seq).replace("-", "").upper()
        if len(s) > len(best):
            best_id, best = f"{f.stem}|{r.id}", s
print(f"EVEscape ref {best_id} L={len(best)}")
if not best:
    raise SystemExit("no AA sequences for EVEscape reference")
Path("$REF_FA").write_text(f">{best_id}\n{best}\n")
PY
  set +e
  "$PYTHON" -u scripts/prepare_evescape.py \
    --csv-path "$EVESCAPE_CSV" \
    --output "$EVESCAPE" \
    --max-seq-len "$MAX_SEQ_LEN" \
    --pathogen other \
    --ref-seq-file "$REF_FA" \
    --min-match-rate 0.55 \
    --no-standardize
  prep_rc=$?
  set -e
  if [[ "$prep_rc" -ne 0 ]]; then
    echo "WARN: prepare_evescape.py failed (rc=$prep_rc) — EVEscape OFF"
    rm -f "$EVESCAPE"
  fi
fi

EV_ARGS=()
if [[ -f "$EVESCAPE" ]]; then
  EV_ARGS+=(--evescape "$EVESCAPE")
  echo "EVEscape: ON ($EVESCAPE)"
  ls -lah "$EVESCAPE"
else
  echo "EVEscape: OFF (no $EVESCAPE). Official Marks summaries have hiv_env_evescape.csv;"
  echo "  tensor build failed or CSV missing. V1–V5 lit_hotspot_mut_frac still ON."
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
    --data "$DATA" \
    --max-seq-len "$MAX_SEQ_LEN" \
    --lit-hotspot-mask "$HIV_MASK" \
    --no-ha-region-metrics \
    --no-spike-domain-metrics \
    --mutation-rate-scale "$MRS" \
    --n-steps "$NSTEPS" \
    --max-trees "$MAX_TREES" \
    --out "$OUT" \
    "${EV_ARGS[@]}"
done
echo "Done: $(date)"
echo "Outputs: checkpoints/eval_enrichment_${CKPT_NAME}_mrs*{${OUT_SUFFIX}}.json"
