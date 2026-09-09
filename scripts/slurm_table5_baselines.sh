#!/bin/bash
#SBATCH --job-name=t5_base
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/table5_baselines_%j.log
#SBATCH --error=logs/table5_baselines_%j.log
#
# Table 5 / C.1 baseline enrichment: pLM + AR + NeutralBD on COVID / H1N1 / HIV.
#
# COVID:
#   VIRUS=covid sbatch --qos=mig-max scripts/slurm_table5_baselines.sh
# H1N1:
#   VIRUS=h1n1 sbatch --qos=mig-max scripts/slurm_table5_baselines.sh
# H3N2:
#   VIRUS=h3n2 METHODS="neutral_bd artreeformer_adapted" \
#     OUT=checkpoints/eval_table5_h3n2_baselines_plmnll.json \
#     sbatch --qos=mig-max scripts/slurm_table5_baselines.sh
# HIV (V1–V5 antigenic mask; no EVEscape tensor):
#   VIRUS=hiv_temporal METHODS="neutral_bd artreeformer_adapted" \
#     sbatch --qos=mig-max scripts/slurm_table5_baselines.sh
#   VIRUS=hiv_geo METHODS="neutral_bd empirical_bd artreeformer_adapted" \
#     sbatch --qos=mig-max scripts/slurm_table5_baselines.sh
# Antibody OAS clones (C.1 free-topo Neutral CTMC / AR pLM NLL):
#   VIRUS=ab_oas METHODS="neutral_bd artreeformer_adapted" \
#     OUT=checkpoints/eval_table5_ab_oas_baselines_plmnll.json \
#     sbatch --qos=mig-max scripts/slurm_table5_baselines.sh
# PLM only (if AR pools missing):
#   METHODS="plm_prior" VIRUS=covid sbatch --qos=mig-max scripts/slurm_table5_baselines.sh

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

VIRUS="${VIRUS:-covid}"
# Leave METHODS unset here so virus-specific defaults apply when not exported.
MAX_TREES="${MAX_TREES:-20}"
N="${N:-16}"

EVESCAPE_ARGS=()
MASK_ARGS=()
if [[ "$VIRUS" == "covid" ]]; then
  DATA=data/covid/test
  TRAIN=data/covid/train
  L=1280
  EVESCAPE=data/covid/evescape_spike_rbd.pt
  MASK=results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt
  OUT="${OUT:-checkpoints/eval_table5_covid_baselines.json}"
  EVESCAPE_ARGS=(--evescape "$EVESCAPE")
  MASK_ARGS=(--lit-hotspot-mask "$MASK")
  METHODS="${METHODS:-plm_prior artreeformer_adapted}"
elif [[ "$VIRUS" == "h1n1" ]]; then
  DATA=data/h1n1/test
  TRAIN=data/h1n1/train
  L=566
  EVESCAPE=data/evescape_h1n1_ha.pt
  MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_h1_nmicrobiol_lit.pt
  OUT="${OUT:-checkpoints/eval_table5_h1n1_baselines.json}"
  EVESCAPE_ARGS=(--evescape "$EVESCAPE")
  MASK_ARGS=(--lit-hotspot-mask "$MASK")
  METHODS="${METHODS:-plm_prior artreeformer_adapted}"
elif [[ "$VIRUS" == "h3n2" ]]; then
  DATA=data/h3n2/test
  TRAIN=data/h3n2/train
  L=566
  MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt
  OUT="${OUT:-checkpoints/eval_table5_h3n2_baselines_plmnll.json}"
  MASK_ARGS=(--lit-hotspot-mask "$MASK")
  METHODS="${METHODS:-neutral_bd artreeformer_adapted}"
elif [[ "$VIRUS" == "hiv_temporal" ]]; then
  DATA=data/hiv_temporal/test
  TRAIN=data/hiv_temporal/train
  L=900
  MASK=results/hiv_env_mask/mut_hotspot_mask_hiv_env_v1v5.pt
  OUT="${OUT:-checkpoints/eval_table5_hiv_temporal_baselines.json}"
  MASK_ARGS=(--lit-hotspot-mask "$MASK")
  METHODS="${METHODS:-neutral_bd artreeformer_adapted}"
elif [[ "$VIRUS" == "hiv_geo" ]]; then
  DATA=data/hiv_geo/test
  TRAIN=data/hiv_geo/train
  L=900
  MASK=results/hiv_env_mask/mut_hotspot_mask_hiv_env_v1v5.pt
  OUT="${OUT:-checkpoints/eval_table5_hiv_geo_baselines.json}"
  MASK_ARGS=(--lit-hotspot-mask "$MASK")
  METHODS="${METHODS:-neutral_bd artreeformer_adapted}"
elif [[ "$VIRUS" == "ab_oas" ]]; then
  DATA=data/ab_clones_1m/test
  TRAIN=data/ab_clones_1m/train
  L=160
  MASK=results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt
  OUT="${OUT:-checkpoints/eval_table5_ab_oas_baselines_plmnll.json}"
  if [[ -f "$MASK" ]]; then
    MASK_ARGS=(--lit-hotspot-mask "$MASK")
  fi
  METHODS="${METHODS:-neutral_bd artreeformer_adapted}"
else
  echo "Unknown VIRUS=$VIRUS (covid|h1n1|h3n2|hiv_temporal|hiv_geo|ab_oas)"; exit 1
fi

echo "Start: $(date) virus=$VIRUS methods=($METHODS) N=$N L=$L out=$OUT"
# shellcheck disable=SC2086
$PYTHON scripts/eval_table5_baselines.py \
  --data "$DATA" --train-data "$TRAIN" \
  --max-seq-len "$L" \
  "${EVESCAPE_ARGS[@]}" \
  "${MASK_ARGS[@]}" \
  --methods $METHODS \
  --N "$N" --max-trees "$MAX_TREES" \
  --out "$OUT"

echo "Done: $(date) -> $OUT"
