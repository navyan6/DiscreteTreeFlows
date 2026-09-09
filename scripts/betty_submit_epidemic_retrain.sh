#!/bin/bash
# Submit epidemic/outbreak v2 retrains (new ckpt dirs only — does not overwrite Wave ckpts).
#
# Prereq on Betty:
#   data/filo_l, data/covid_epidemic, data/h3n2_epidemic, data/h1n1_epidemic
#   bash scripts/betty_prep_epidemic_splits.sh
#
# Usage:
#   bash scripts/betty_submit_epidemic_retrain.sh          # COVID + flu epidemic (default)
#   bash scripts/betty_submit_epidemic_retrain.sh --with-filo
#   bash scripts/betty_submit_epidemic_retrain.sh --filo-only

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs benchmarks/results/tables

FILO=0
COVID=1
FLU=1
TL_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --filo-only) COVID=0; FLU=0; FILO=1 ;;
    --with-filo) FILO=1 ;;
    --skip-filo) FILO=0 ;;
    --tl-only) COVID=0; FLU=0; TL_ONLY=1; FILO=1 ;;
  esac
done

IDS_JSON=benchmarks/results/tables/epidemic_retrain_job_ids.json
echo "{" > "$IDS_JSON.tmp"
first=1

submit_chain() {
  local name="$1"
  local data_root="$2"
  local max_len="$3"
  local ckpt="$4"
  local hotspot="${5:-0}"
  local mask="${6:-}"
  local init_ckpt="${7:-}"
  local lr="${8:-}"

  local j_pipe j_pre j_train
  j_pipe=$(sbatch --parsable \
    --export=ALL,DATA_ROOT="$data_root" \
    scripts/slurm_epidemic_pipeline.sh)
  echo "  $name pipeline -> $j_pipe"

  j_pre=$(sbatch --parsable --dependency=afterok:"$j_pipe" \
    --export=ALL,DATA_ROOT="$data_root",MAX_SEQ_LEN="$max_len" \
    scripts/slurm_epidemic_precompute.sh)
  echo "  $name precompute -> $j_pre"

  local extra="ALL,DATA_ROOT=$data_root,CKPT_DIR=$ckpt,MAX_SEQ_LEN=$max_len,HOTSPOT=$hotspot"
  if [[ -n "$mask" ]]; then
    extra="$extra,MUT_HOTSPOT_MASK=$mask"
  fi
  if [[ -n "$init_ckpt" ]]; then
    extra="$extra,INIT_CHECKPOINT=$init_ckpt"
    if [[ -n "$lr" ]]; then
      extra="$extra,LR=$lr"
    fi
  fi
  j_train=$(sbatch --parsable --dependency=afterok:"$j_pre" --qos=mig-max \
    --job-name="ep_${name}" --export="$extra" \
    scripts/slurm_epidemic_train_mut_recovery.sh)
  echo "  $name train -> $j_train"

  if [[ $first -eq 1 ]]; then first=0; else echo "," >> "$IDS_JSON.tmp"; fi
  printf '  "%s_pipeline": %s,\n  "%s_precompute": %s,\n  "%s_train": %s' \
    "$name" "$j_pipe" "$name" "$j_pre" "$name" "$j_train" >> "$IDS_JSON.tmp"
}

comma() {
  if [[ $first -eq 0 ]]; then echo "," >> "$IDS_JSON.tmp"; fi
}

if [[ "$FILO" == "1" ]]; then
  if [[ -d data/filo_l/train ]] && ls data/filo_l/train/*_group_*.fasta &>/dev/null; then
    echo "=== Filo L outbreak ==="
    if [[ "$TL_ONLY" != "1" ]]; then
      submit_chain "filo_l_v1" "data/filo_l" "900" "checkpoints/filo_l_v1_mutrec" "0" ""
      comma
      submit_chain "filo_l_v1_lit" "data/filo_l" "900" "checkpoints/filo_l_v1_lit_mutrec" "1" \
        "results/bdbv_l_mask/mut_hotspot_mask_lit.pt"
    fi
    TL_INIT=""
    for cand in checkpoints/covid_v5_mutrec/best.pt checkpoints/h3n2_v3_lit_hotspot/best.pt; do
      if [[ -f "$cand" ]]; then TL_INIT="$cand"; break; fi
    done
    if [[ -n "$TL_INIT" ]]; then
      echo "=== Filo L TL (init $TL_INIT) ==="
      comma
      submit_chain "filo_l_v1_tl" "data/filo_l" "900" "checkpoints/filo_l_v1_tl_mutrec" "0" "" \
        "$TL_INIT" "3e-5"
    else
      echo "SKIP filo_l_v1_tl (no covid/flu init checkpoint on disk)"
    fi
  else
    echo "SKIP filo_l (no groups)"
  fi
fi

if [[ "$COVID" == "1" ]]; then
  if [[ -d data/covid_epidemic/train ]] && ls data/covid_epidemic/train/*_group_*.fasta &>/dev/null; then
    echo "=== COVID epidemic ==="
    comma
    submit_chain "covid_v6_epidemic" "data/covid_epidemic" "1280" \
      "checkpoints/covid_v6_epidemic_mutrec" "0" ""
  fi
fi

if [[ "$FLU" == "1" ]]; then
  if [[ -d data/h3n2_epidemic/train ]] && ls data/h3n2_epidemic/train/*_group_*.fasta &>/dev/null; then
    echo "=== H3N2 epidemic ==="
    comma
    submit_chain "h3n2_v4_epidemic" "data/h3n2_epidemic" "566" \
      "checkpoints/h3n2_v4_epidemic_mutrec" "0" ""
    comma
    submit_chain "h3n2_v4_epidemic_lit" "data/h3n2_epidemic" "566" \
      "checkpoints/h3n2_v4_epidemic_lit_mutrec" "1" \
      "results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt"
  fi
  if [[ -d data/h1n1_epidemic/train ]] && ls data/h1n1_epidemic/train/*_group_*.fasta &>/dev/null; then
    echo "=== H1N1 epidemic ==="
    comma
    submit_chain "h1n1_v3_epidemic" "data/h1n1_epidemic" "566" \
      "checkpoints/h1n1_v3_epidemic_mutrec" "0" ""
    comma
    submit_chain "h1n1_v3_epidemic_lit" "data/h1n1_epidemic" "566" \
      "checkpoints/h1n1_v3_epidemic_lit_mutrec" "1" \
      "results/flu_mutfreq_vs_lit/mut_hotspot_mask_h1_nmicrobiol_lit.pt"
  fi
fi

echo "" >> "$IDS_JSON.tmp"
echo "}" >> "$IDS_JSON.tmp"
mv "$IDS_JSON.tmp" "$IDS_JSON"
echo "Wrote $IDS_JSON"
cat "$IDS_JSON"
