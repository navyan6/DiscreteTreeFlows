#!/bin/bash
#SBATCH --job-name=h3n2_lithot
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/h3n2_lithot_%j.log
#SBATCH --error=logs/h3n2_lithot_%j.log
#
# H3N2 lit-hotspot train (nmicrobiol antigenic sites).
# STRICT: MUT_HOTSPOT_MASK must be under results/flu_mutfreq_vs_lit/ (never covid/PMC).
#
# Example:
#   HOTSPOT=1 CKPT_DIR=checkpoints/h3n2_v3_lit_hotspot \
#     MUT_HOTSPOT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt \
#     MUT_HOTSPOT_FRAC= MUT_HOTSPOT_TOPK= \
#     MUT_HOTSPOT_WEIGHT=5 MUT_HOTSPOT_FORCE=1 \
#     sbatch --qos=mig-max --export=ALL,HOTSPOT,CKPT_DIR,MUT_HOTSPOT_MASK,MUT_HOTSPOT_FRAC,MUT_HOTSPOT_TOPK,MUT_HOTSPOT_WEIGHT,MUT_HOTSPOT_FORCE \
#       scripts/slurm_h3n2_train_lit_hotspot.sh

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}

cd ~/DiscreteTreeFlows

CKPT_DIR="${CKPT_DIR:-checkpoints/h3n2_v3_lit_hotspot}"
LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
MUT_NORMALIZE="${MUT_NORMALIZE:-count}"
ENTROPY_ALPHA="${ENTROPY_ALPHA:-3.0}"
ENTROPY_ALPHA_CONS="${ENTROPY_ALPHA_CONS:-1.0}"
PSSM_GATE="${PSSM_GATE:-1}"
MUT_AA_EMB="${MUT_AA_EMB:-1}"
MUT_AA_EMB_DIM="${MUT_AA_EMB_DIM:-16}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"

MUT_HOTSPOT_MASK="${MUT_HOTSPOT_MASK:-results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt}"
MUT_HOTSPOT_FRAC="${MUT_HOTSPOT_FRAC:-}"
MUT_HOTSPOT_TOPK="${MUT_HOTSPOT_TOPK:-}"
MUT_HOTSPOT_WEIGHT="${MUT_HOTSPOT_WEIGHT:-5.0}"
MUT_HOTSPOT_FORCE="${MUT_HOTSPOT_FORCE:-1}"
MUT_HOTSPOT_SCORE="${MUT_HOTSPOT_SCORE:-entropy}"
FITNESS_BETA="${FITNESS_BETA:-0}"
FITNESS_SCORE="${FITNESS_SCORE:-log_R0}"

# Refuse covid/PMC masks — H3N2 lit only.
case "$MUT_HOTSPOT_MASK" in
  *covid_mutfreq*|*pmc_lit*|*/covid/*)
    echo "REFUSING covid/PMC hotspot mask for H3N2 train: $MUT_HOTSPOT_MASK" >&2
    exit 2
    ;;
esac
if [[ ! -f "$MUT_HOTSPOT_MASK" ]]; then
  echo "MISSING flu hotspot mask: $MUT_HOTSPOT_MASK" >&2
  exit 1
fi

mkdir -p logs "$CKPT_DIR"

EXTRA_ARGS=()
if [[ "$PSSM_GATE" == "1" || "$PSSM_GATE" == "true" ]]; then
    EXTRA_ARGS+=(--pssm-gate)
fi
if [[ "$MUT_AA_EMB" == "1" || "$MUT_AA_EMB" == "true" ]]; then
    EXTRA_ARGS+=(--mut-aa-emb --mut-aa-emb-dim "$MUT_AA_EMB_DIM")
fi
if [[ -n "$ENTROPY_ALPHA_CONS" ]]; then
    EXTRA_ARGS+=(--entropy-weight-alpha-cons "$ENTROPY_ALPHA_CONS")
fi
if [[ "$FITNESS_BETA" != "0" && "$FITNESS_BETA" != "0.0" ]]; then
    EXTRA_ARGS+=(--fitness-beta "$FITNESS_BETA" --fitness-score "$FITNESS_SCORE")
fi

# Prefer explicit flu lit mask; clear FRAC/TOPK when mask is set.
if [[ -n "$MUT_HOTSPOT_MASK" ]]; then
    EXTRA_ARGS+=(--mut-hotspot-mask "$MUT_HOTSPOT_MASK")
    EXTRA_ARGS+=(--mut-hotspot-weight "$MUT_HOTSPOT_WEIGHT")
elif [[ -n "${MUT_HOTSPOT_TOPK:-}" ]]; then
    EXTRA_ARGS+=(--mut-hotspot-topk "$MUT_HOTSPOT_TOPK")
    EXTRA_ARGS+=(--mut-hotspot-weight "$MUT_HOTSPOT_WEIGHT")
    EXTRA_ARGS+=(--mut-hotspot-score "$MUT_HOTSPOT_SCORE")
elif [[ -n "${MUT_HOTSPOT_FRAC:-}" ]]; then
    EXTRA_ARGS+=(--mut-hotspot-frac "$MUT_HOTSPOT_FRAC")
    EXTRA_ARGS+=(--mut-hotspot-weight "$MUT_HOTSPOT_WEIGHT")
    EXTRA_ARGS+=(--mut-hotspot-score "$MUT_HOTSPOT_SCORE")
fi
if [[ "${MUT_HOTSPOT_FORCE:-0}" == "1" || "${MUT_HOTSPOT_FORCE:-0}" == "true" ]]; then
    EXTRA_ARGS+=(--mut-hotspot-force)
fi

echo "Start: $(date)"
echo "ckpt_dir=$CKPT_DIR  max_seq_len=$MAX_SEQ_LEN"
echo "data=data/h3n2/{train,val,test}"
echo "lambda_mut=$LAMBDA_MUT  lambda_cons=$LAMBDA_CONS  mut_normalize=$MUT_NORMALIZE"
echo "lambda_semi=$LAMBDA_SEMI  pssm_gate=$PSSM_GATE  mut_aa_emb=$MUT_AA_EMB"
echo "entropy_alpha=$ENTROPY_ALPHA  entropy_alpha_cons=${ENTROPY_ALPHA_CONS:-same}"
echo "fitness_beta=$FITNESS_BETA  fitness_score=$FITNESS_SCORE"
echo "mut_hotspot_mask=$MUT_HOTSPOT_MASK"
echo "mut_hotspot_topk=${MUT_HOTSPOT_TOPK:-}  mut_hotspot_frac=${MUT_HOTSPOT_FRAC:-}"
echo "mut_hotspot_weight=$MUT_HOTSPOT_WEIGHT  mut_hotspot_force=$MUT_HOTSPOT_FORCE"
echo "site_levers: --per-site-pos-emb --use-site-entropy --use-entropy-loss-weighting --use-entropy-cons-weighting (hardcoded ON); pssm_gate=$PSSM_GATE mut_aa_emb=$MUT_AA_EMB"
ls -la "$MUT_HOTSPOT_MASK"
echo "Loading flu lit hotspot mask: $MUT_HOTSPOT_MASK"

$PYTHON -u scripts/train.py \
    --data        data/h3n2/train \
    --val-data    data/h3n2/val \
    --test-data   data/h3n2/test \
    --max-seq-len "$MAX_SEQ_LEN" \
    --epochs      100 \
    --patience    50 \
    --bridge-c    1.0 \
    --lambda-mut  "$LAMBDA_MUT" \
    --lambda-cons "$LAMBDA_CONS" \
    --mut-normalize "$MUT_NORMALIZE" \
    --lambda-semi "$LAMBDA_SEMI" \
    --per-site-pos-emb \
    --use-site-entropy \
    --use-entropy-loss-weighting \
    --use-entropy-cons-weighting \
    --entropy-source empirical \
    --entropy-weight-alpha "$ENTROPY_ALPHA" \
    --entropy-weight-floor 1.0 \
    --ckpt-dir    "$CKPT_DIR" \
    --resume \
    "${EXTRA_ARGS[@]}"

echo "Done: $(date)"
