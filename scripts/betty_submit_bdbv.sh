#!/bin/bash
# End-to-end BDBV L TreeSBM: ingest → window → splits → SLURM pipeline/train/eval.
#
# Run on Betty login node:
#   bash scripts/betty_submit_bdbv.sh                 # SLURM only (default; data must exist)
#   bash scripts/betty_submit_bdbv.sh --ingest        # NCBI → splits → SLURM
#   bash scripts/betty_submit_bdbv.sh --ingest-only   # data pipeline only (login node)
#   bash scripts/betty_submit_bdbv.sh --eval-only     # after Pathoplexus merge or re-eval
#
# Preflight (no SSH required to edit; run on Betty before submit):
#   bash scripts/bdbv_preflight.sh
#
# Never overwrites paper ckpts (covid_v5_mutrec, h1n1_v2, h3n2_v3, hiv_geo_v1).

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
PYTHON="${PYTHON:-/vast/home/n/nnori/.conda/envs/treesbm/bin/python}"

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs benchmarks/results/tables data/bdbv/references

INGEST=0
SUBMIT=1
EVAL_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --ingest) INGEST=1 ;;
    --ingest-only) SUBMIT=0; INGEST=1 ;;
    --eval-only) INGEST=0; EVAL_ONLY=1 ;;
    --no-submit) SUBMIT=0 ;;
  esac
done

if [[ "$SUBMIT" == "1" && "$EVAL_ONLY" != "1" ]]; then
  bash scripts/bdbv_preflight.sh
fi

IDS_JSON=benchmarks/results/tables/bdbv_job_ids.json
echo "{" > "$IDS_JSON.tmp"
first=1

submit() {
  local name="$1"; shift
  local jid
  jid=$(sbatch --parsable "$@")
  echo "submitted $name -> $jid"
  if [[ $first -eq 1 ]]; then first=0; else echo "," >> "$IDS_JSON.tmp"; fi
  printf '  "%s": %s' "$name" "$jid" >> "$IDS_JSON.tmp"
  REPLY="$jid"
}

if [[ "$INGEST" == "1" ]]; then
  echo "=== Phase 0: NCBI + Nextstrain ==="
  $PYTHON scripts/download_ebolavirus_ncbi.py --max-per-species 2000 || true
  $PYTHON scripts/download_ebolavirus_nextstrain.py || true

  echo "=== Phase 1: L extract + window + lit mask ==="
  $PYTHON scripts/bdbv_extract_l.py
  $PYTHON scripts/bdbv_define_l_window.py
  $PYTHON scripts/build_bdbv_l_lit_hotspot_mask.py

  echo "=== Phase 2: splits (outbreak primary) ==="
  $PYTHON scripts/prepare_filo_outbreak.py --out-base data/filo_l
  $PYTHON scripts/prepare_filo_outbreak.py --track-b --test-outbreak ebov_wa_2013_2016 \
    --out-base data/filo_l_track_b || true
  $PYTHON scripts/audit_filo_splits.py --data data/filo_l
  # Legacy temporal ablation (optional comparison)
  $PYTHON scripts/prepare_bdbv_temporal.py || true
  $PYTHON scripts/prepare_bdbv_temporal.py --pan-ebolavirus --out-base data/bdbv_pan_temporal || true
  $PYTHON scripts/audit_bdbv_splits.py || true
fi

if [[ "$SUBMIT" != "1" ]]; then
  echo "" >> "$IDS_JSON.tmp"
  echo "}" >> "$IDS_JSON.tmp"
  mv "$IDS_JSON.tmp" "$IDS_JSON"
  echo "Ingest done (no SLURM jobs submitted)."
  exit 0
fi

if [[ "$EVAL_ONLY" == "1" ]]; then
  for CKPT in filo_l_v1_mutrec filo_l_v1_lit_mutrec bdbv_v1_mutrec bdbv_pan_v1_mutrec; do
    P="checkpoints/${CKPT}/best.pt"
    [[ -f "$P" ]] || continue
    DR=data/filo_l
    [[ "$CKPT" == bdbv_* ]] && DR=data/bdbv_temporal
    [[ "$CKPT" == bdbv_pan_* ]] && DR=data/bdbv_pan_temporal
    submit "eval_${CKPT}" --qos=mig-max --job-name="ev_${CKPT}" \
      --export=ALL,CKPT="$P",DATA_ROOT="$DR" \
      scripts/slurm_eval_bdbv.sh
    submit "cov_${CKPT}" --qos=mig-max --job-name="cov_${CKPT}" \
      --export=ALL,CHECKPOINT="$P",DATA_ROOT="$DR" \
      scripts/slurm_bdbv_coverage.sh
  done
  echo "" >> "$IDS_JSON.tmp"
  echo "}" >> "$IDS_JSON.tmp"
  mv "$IDS_JSON.tmp" "$IDS_JSON"
  echo "Eval-only jobs -> $IDS_JSON"
  exit 0
fi

echo "=== Phase 3: tree pipeline (filo_l outbreak — primary) ==="
J_PIPE=$(sbatch --parsable \
  --export=ALL,DATA_ROOT=data/filo_l \
  scripts/slurm_bdbv_pipeline.sh)
echo "pipeline filo_l -> $J_PIPE"
first=0
printf '  "filo_l_pipeline": %s' "$J_PIPE" >> "$IDS_JSON.tmp"

J_PRE=$(sbatch --parsable --dependency=afterok:"$J_PIPE" \
  --export=ALL,DATA_ROOT=data/filo_l \
  scripts/slurm_bdbv_precompute.sh)
echo "precompute filo_l -> $J_PRE"
echo "," >> "$IDS_JSON.tmp"
printf '  "filo_l_precompute": %s' "$J_PRE" >> "$IDS_JSON.tmp"

DEP="afterok:$J_PRE"

submit "filo_l_v1_mutrec" --dependency="$DEP" --qos=mig-max --job-name=filo_v1 \
  --export=ALL,CKPT_DIR=checkpoints/filo_l_v1_mutrec,DATA_ROOT=data/filo_l,HOTSPOT=0 \
  scripts/slurm_bdbv_train_mut_recovery.sh
J_FILO_V1="$REPLY"

submit "filo_l_v1_lit_mutrec" --dependency="$DEP" --qos=mig-max --job-name=filo_v1_lit \
  --export=ALL,CKPT_DIR=checkpoints/filo_l_v1_lit_mutrec,DATA_ROOT=data/filo_l,HOTSPOT=1 \
  scripts/slurm_bdbv_train_mut_recovery.sh
J_FILO_V1_LIT="$REPLY"

submit "eval_filo_l_v1" --dependency="afterok:$J_FILO_V1" --qos=mig-max --job-name=ev_filo_v1 \
  --export=ALL,CKPT=checkpoints/filo_l_v1_mutrec/best.pt,DATA_ROOT=data/filo_l \
  scripts/slurm_eval_bdbv.sh

submit "eval_filo_l_v1_lit" --dependency="afterok:$J_FILO_V1_LIT" --qos=mig-max --job-name=ev_filo_lit \
  --export=ALL,CKPT=checkpoints/filo_l_v1_lit_mutrec/best.pt,DATA_ROOT=data/filo_l \
  scripts/slurm_eval_bdbv.sh

submit "cov_filo_l_v1" --dependency="afterok:$J_FILO_V1" --qos=mig-max --job-name=cov_filo_v1 \
  --export=ALL,CHECKPOINT=checkpoints/filo_l_v1_mutrec/best.pt,DATA_ROOT=data/filo_l \
  scripts/slurm_bdbv_coverage.sh

submit "ftopo_filo_l_v1" --dependency="afterok:$J_FILO_V1" --qos=mig-max --job-name=ft_filo \
  --export=ALL,VIRUS=bdbv,CKPT=checkpoints/filo_l_v1_mutrec/best.pt \
  scripts/slurm_eval_fixed_topo.sh

# Legacy temporal ablation (optional; skip if ingest not run)
if [[ -d data/bdbv_temporal/train ]] && ls data/bdbv_temporal/train/*_group_*.fasta &>/dev/null; then
  echo "=== Ablation: bdbv temporal ==="
  J_PIPE_B=$(sbatch --parsable scripts/slurm_bdbv_pipeline.sh)
  echo "," >> "$IDS_JSON.tmp"
  printf '  "bdbv_pipeline": %s' "$J_PIPE_B" >> "$IDS_JSON.tmp"
  J_PRE_B=$(sbatch --parsable --dependency=afterok:"$J_PIPE_B" \
    --export=ALL,DATA_ROOT=data/bdbv_temporal \
    scripts/slurm_bdbv_precompute.sh)
  echo "," >> "$IDS_JSON.tmp"
  printf '  "bdbv_precompute": %s' "$J_PRE_B" >> "$IDS_JSON.tmp"
  DEP_B="afterok:$J_PRE_B"
  submit "bdbv_v1_mutrec" --dependency="$DEP_B" --qos=mig-max --job-name=bdbv_v1 \
    --export=ALL,CKPT_DIR=checkpoints/bdbv_v1_mutrec,DATA_ROOT=data/bdbv_temporal,HOTSPOT=0 \
    scripts/slurm_bdbv_train_mut_recovery.sh || true
fi

# Pan temporal ablation
if [[ -d data/bdbv_pan_temporal/train ]] && ls data/bdbv_pan_temporal/train/*_group_*.fasta &>/dev/null; then
  J_PIPE_PAN=$(sbatch --parsable \
    --export=ALL,DATA_ROOT=data/bdbv_pan_temporal \
    scripts/slurm_bdbv_pipeline.sh)
  echo "," >> "$IDS_JSON.tmp"
  printf '  "bdbv_pan_pipeline": %s' "$J_PIPE_PAN" >> "$IDS_JSON.tmp"
  J_PRE_PAN=$(sbatch --parsable --dependency=afterok:"$J_PIPE_PAN" \
    --export=ALL,DATA_ROOT=data/bdbv_pan_temporal \
    scripts/slurm_bdbv_precompute.sh)
  echo "," >> "$IDS_JSON.tmp"
  printf '  "bdbv_pan_precompute": %s' "$J_PRE_PAN" >> "$IDS_JSON.tmp"
  DEP_PAN="afterok:$J_PRE_PAN"
  submit "bdbv_pan_v1_mutrec" --dependency="$DEP_PAN" --qos=mig-max --job-name=bdbv_pan_v1 \
    --export=ALL,CKPT_DIR=checkpoints/bdbv_pan_v1_mutrec,DATA_ROOT=data/bdbv_pan_temporal,HOTSPOT=0 \
    scripts/slurm_bdbv_train_mut_recovery.sh || true
fi

echo "" >> "$IDS_JSON.tmp"
echo "}" >> "$IDS_JSON.tmp"
mv "$IDS_JSON.tmp" "$IDS_JSON"
echo "Wrote $IDS_JSON"
cat "$IDS_JSON"
exit 0
