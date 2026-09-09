#!/bin/bash
#SBATCH --job-name=covid_mutrec
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/covid_mutrec_%j.log
#SBATCH --error=logs/covid_mutrec_%j.log
#
# SH26 fast-path train (covid_v5): stronger mut loss + PSSM gate + mut AA emb.
#   --per-site-pos-emb --lambda-semi
#   --lambda-mut 12 --lambda-cons 0.5 --mut-normalize count
#   --entropy-weight-alpha-cons 1.0 (mut alpha stays 3.0)
#   --pssm-gate --mut-aa-emb
#
# covid_v6_hotspot (HOTSPOT=1): v5 recipe + hard top-frac MSA entropy hotspots
#   --mut-hotspot-frac 0.15 --mut-hotspot-weight 5.0
#   (equiv. --mut-hotspot-topk 192 for L=1280; use topk=80 for a tighter set)
#
# Default CKPT_DIR=checkpoints/covid_v5_mutrec (new arch → do NOT resume v4).
# Set V4_COMPAT=1 to restore the older v4 recipe into covid_v4_mutrec.
# Set HOTSPOT=1 for covid_v6_hotspot (do NOT resume into v5).
#
# Env overrides:
#   sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh
#   LAMBDA_MUT=16 LAMBDA_CONS=0.25 sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh
#   DEEP_MUT=1 CKPT_DIR=checkpoints/covid_v5_deepmut \
#     sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh
#   V4_COMPAT=1 CKPT_DIR=checkpoints/covid_v4_mutrec \
#     sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh
#   HOTSPOT=1 CKPT_DIR=checkpoints/covid_v6_hotspot \
#     sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh
#   HOTSPOT=1 MUT_HOTSPOT_TOPK=80 MUT_HOTSPOT_FRAC= \
#     CKPT_DIR=checkpoints/covid_v6_hotspot_topk80 \
#     sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh
#   # Preferred antigenic train (curated PMC lit mask + force):
#   HOTSPOT=1 CKPT_DIR=checkpoints/covid_v7_pmc_hotspot \
#     MUT_HOTSPOT_MASK=results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt \
#     MUT_HOTSPOT_FRAC= MUT_HOTSPOT_TOPK= \
#     MUT_HOTSPOT_WEIGHT=5 MUT_HOTSPOT_FORCE=1 \
#     sbatch --qos=mig-max --export=ALL,HOTSPOT,CKPT_DIR,MUT_HOTSPOT_MASK,MUT_HOTSPOT_FRAC,MUT_HOTSPOT_TOPK,MUT_HOTSPOT_WEIGHT,MUT_HOTSPOT_FORCE \
#       scripts/slurm_covid_train_mut_recovery.sh
#   FITNESS_BETA=1.0 CKPT_DIR=checkpoints/covid_v5_fitness \
#     sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh
#     # §4.2 Option A R0 tilt; also FITNESS_SCORE=log_softmax

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}

cd ~/DiscreteTreeFlows

V4_COMPAT="${V4_COMPAT:-0}"
HOTSPOT="${HOTSPOT:-0}"
if [[ "$V4_COMPAT" == "1" || "$V4_COMPAT" == "true" ]]; then
    CKPT_DIR="${CKPT_DIR:-checkpoints/covid_v4_mutrec}"
    LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
    LAMBDA_MUT="${LAMBDA_MUT:-8.0}"
    LAMBDA_CONS="${LAMBDA_CONS:-1.0}"
    MUT_NORMALIZE="${MUT_NORMALIZE:-mean}"
    ENTROPY_ALPHA="${ENTROPY_ALPHA:-3.0}"
    ENTROPY_ALPHA_CONS="${ENTROPY_ALPHA_CONS:-}"
    PSSM_GATE="${PSSM_GATE:-0}"
    MUT_AA_EMB="${MUT_AA_EMB:-0}"
elif [[ "$HOTSPOT" == "1" || "$HOTSPOT" == "true" ]]; then
    CKPT_DIR="${CKPT_DIR:-checkpoints/covid_v6_hotspot}"
    LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
    LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
    LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
    MUT_NORMALIZE="${MUT_NORMALIZE:-count}"
    ENTROPY_ALPHA="${ENTROPY_ALPHA:-3.0}"
    ENTROPY_ALPHA_CONS="${ENTROPY_ALPHA_CONS:-1.0}"
    PSSM_GATE="${PSSM_GATE:-1}"
    MUT_AA_EMB="${MUT_AA_EMB:-1}"
    # Default: top 15% of spike columns (~192 of 1280). Override with TOPK
    # or MUT_HOTSPOT_MASK (curated PMC lit). Empty FRAC/TOPK when mask is set.
    MUT_HOTSPOT_MASK="${MUT_HOTSPOT_MASK:-}"
    if [[ -n "$MUT_HOTSPOT_MASK" ]]; then
        MUT_HOTSPOT_FRAC="${MUT_HOTSPOT_FRAC:-}"
        MUT_HOTSPOT_TOPK="${MUT_HOTSPOT_TOPK:-}"
    else
        MUT_HOTSPOT_FRAC="${MUT_HOTSPOT_FRAC:-0.15}"
        MUT_HOTSPOT_TOPK="${MUT_HOTSPOT_TOPK:-}"
    fi
    MUT_HOTSPOT_WEIGHT="${MUT_HOTSPOT_WEIGHT:-5.0}"
    MUT_HOTSPOT_FORCE="${MUT_HOTSPOT_FORCE:-0}"
else
    CKPT_DIR="${CKPT_DIR:-checkpoints/covid_v5_mutrec}"
    LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
    LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
    LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
    MUT_NORMALIZE="${MUT_NORMALIZE:-count}"
    ENTROPY_ALPHA="${ENTROPY_ALPHA:-3.0}"
    ENTROPY_ALPHA_CONS="${ENTROPY_ALPHA_CONS:-1.0}"
    PSSM_GATE="${PSSM_GATE:-1}"
    MUT_AA_EMB="${MUT_AA_EMB:-1}"
    MUT_HOTSPOT_MASK="${MUT_HOTSPOT_MASK:-}"
    MUT_HOTSPOT_FRAC="${MUT_HOTSPOT_FRAC:-}"
    MUT_HOTSPOT_TOPK="${MUT_HOTSPOT_TOPK:-}"
    MUT_HOTSPOT_WEIGHT="${MUT_HOTSPOT_WEIGHT:-5.0}"
    MUT_HOTSPOT_FORCE="${MUT_HOTSPOT_FORCE:-0}"
fi

mkdir -p logs "$CKPT_DIR"

DEEP_MUT="${DEEP_MUT:-0}"
MUT_AA_EMB_DIM="${MUT_AA_EMB_DIM:-16}"
FITNESS_BETA="${FITNESS_BETA:-0}"
FITNESS_SCORE="${FITNESS_SCORE:-log_R0}"

EXTRA_ARGS=()
if [[ "$DEEP_MUT" == "1" || "$DEEP_MUT" == "true" ]]; then
    EXTRA_ARGS+=(--deep-mut-head)
fi
if [[ "$PSSM_GATE" == "1" || "$PSSM_GATE" == "true" ]]; then
    EXTRA_ARGS+=(--pssm-gate)
fi
if [[ "$MUT_AA_EMB" == "1" || "$MUT_AA_EMB" == "true" ]]; then
    EXTRA_ARGS+=(--mut-aa-emb --mut-aa-emb-dim "$MUT_AA_EMB_DIM")
fi
if [[ -n "$ENTROPY_ALPHA_CONS" ]]; then
    EXTRA_ARGS+=(--entropy-weight-alpha-cons "$ENTROPY_ALPHA_CONS")
fi
# §4.2 Option A fitness tilt (0 = off / backward compatible)
if [[ "$FITNESS_BETA" != "0" && "$FITNESS_BETA" != "0.0" ]]; then
    EXTRA_ARGS+=(--fitness-beta "$FITNESS_BETA" --fitness-score "$FITNESS_SCORE")
fi

# Hard MSA hotspots (opt-in; default on under HOTSPOT=1).
# MUT_HOTSPOT_SCORE=entropy (default) | mut_freq  (MSA-select → tree-apply)
# MUT_HOTSPOT_MASK=path.pt overrides score ranking (preferred for PMC lit).
MUT_HOTSPOT_SCORE="${MUT_HOTSPOT_SCORE:-entropy}"
MUT_HOTSPOT_MASK="${MUT_HOTSPOT_MASK:-}"
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
echo "ckpt_dir=$CKPT_DIR  v4_compat=$V4_COMPAT  hotspot=$HOTSPOT"
echo "lambda_mut=$LAMBDA_MUT  lambda_cons=$LAMBDA_CONS  mut_normalize=$MUT_NORMALIZE"
echo "lambda_semi=$LAMBDA_SEMI  pssm_gate=$PSSM_GATE  mut_aa_emb=$MUT_AA_EMB  deep_mut=$DEEP_MUT"
echo "entropy_alpha=$ENTROPY_ALPHA  entropy_alpha_cons=${ENTROPY_ALPHA_CONS:-same}"
echo "fitness_beta=$FITNESS_BETA  fitness_score=$FITNESS_SCORE"
echo "mut_hotspot_mask=${MUT_HOTSPOT_MASK:-}  mut_hotspot_topk=${MUT_HOTSPOT_TOPK:-}  mut_hotspot_frac=${MUT_HOTSPOT_FRAC:-}"
echo "mut_hotspot_score=${MUT_HOTSPOT_SCORE:-}  mut_hotspot_weight=${MUT_HOTSPOT_WEIGHT:-}  "
echo "mut_hotspot_force=${MUT_HOTSPOT_FORCE:-0}"
echo "lambda_br=${LAMBDA_BR:-0.1}  (Table 8: set LAMBDA_BR=0 to drop BL head loss)"
echo "site_levers: --per-site-pos-emb --use-site-entropy --use-entropy-loss-weighting --use-entropy-cons-weighting (hardcoded ON)"

LAMBDA_BR="${LAMBDA_BR:-0.1}"
LAMBDA_TOP="${LAMBDA_TOP:-0.1}"

$PYTHON -u scripts/train.py \
    --data        data/covid/train \
    --val-data    data/covid/val \
    --test-data   data/covid/test \
    --max-seq-len 1280 \
    --epochs      100 \
    --patience    50 \
    --bridge-c    1.0 \
    --lambda-mut  "$LAMBDA_MUT" \
    --lambda-cons "$LAMBDA_CONS" \
    --lambda-br   "$LAMBDA_BR" \
    --lambda-top  "$LAMBDA_TOP" \
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
