#!/bin/bash
#SBATCH --job-name=eval_eve
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/eval_eve_%j.log
#SBATCH --error=logs/eval_eve_%j.log
#
# EVE (Marks / OATML) evolutionary-index baseline vs TreeSBM generations.
# Metrics (eval_eve_baseline.py): mut/cons recovery + mean EVE of recovered GT
# mutations vs random AA @ mut sites + Pearson/Spearman of TreeSBM mut log-probs
# vs EVE at GT mutating sites.
#
# Requires a precomputed --eve-scores .pt ([L,20], AA order ACDEFGHIKLMNPQRSTVWY).
# Building that tensor is offline (prepare_eve_scores.py + EVE repo); scoring
# itself is cheap — GPU here is for TreeSBM + ESM tree generation.
# Override QOS: sbatch --qos=<qos> scripts/slurm_eval_eve.sh ...
#
# Usage:
#   sbatch scripts/slurm_eval_eve.sh <checkpoint> <data_dir> <eve_scores.pt> [max_seq_len] [out_json]
#
# Examples:
#   # COVID Spike
#   sbatch scripts/slurm_eval_eve.sh \
#     checkpoints/covid_v3_cons/best.pt data/covid/test data/covid/eve_spike.pt 1280
#
#   # H3N2 HA
#   sbatch scripts/slurm_eval_eve.sh \
#     checkpoints/h3n2_v2/best.pt data/h3n2/test data/h3n2/eve_ha.pt 566

set -e

CHECKPOINT="${1:?usage: sbatch slurm_eval_eve.sh <checkpoint> <data_dir> <eve_scores.pt> [max_seq_len] [out_json]}"
DATA_DIR="${2:?usage: sbatch slurm_eval_eve.sh <checkpoint> <data_dir> <eve_scores.pt> [max_seq_len] [out_json]}"
EVE_SCORES="${3:?usage: sbatch slurm_eval_eve.sh <checkpoint> <data_dir> <eve_scores.pt> [max_seq_len] [out_json]}"
MAX_SEQ_LEN="${4:-566}"
OUT_JSON="${5:-}"

mkdir -p logs checkpoints

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"
export TORCH_HOME="${TORCH_HOME:-$LABHOME/torch_cache}"

cd ~/DiscreteTreeFlows

# Prefer lab backup if local ckpt missing
if [[ ! -f "$CHECKPOINT" && -f "$LABHOME/checkpoints_backup/$(basename "$(dirname "$CHECKPOINT")")/best.pt" ]]; then
    name="$(basename "$(dirname "$CHECKPOINT")")"
    mkdir -p "checkpoints/${name}"
    cp -n "$LABHOME/checkpoints_backup/${name}/best.pt" "checkpoints/${name}/best.pt"
    CHECKPOINT="checkpoints/${name}/best.pt"
fi

if [[ ! -f "$CHECKPOINT" ]]; then
    echo "ERROR: checkpoint not found: $CHECKPOINT"
    exit 1
fi
if [[ ! -f "$EVE_SCORES" ]]; then
    echo "ERROR: EVE scores missing: $EVE_SCORES"
    echo "Build offline (see EXTERNAL.md § EVE / scripts/prepare_eve_scores.py):"
    echo "  git clone https://github.com/OATML-Markslab/EVE \$LABHOME/baselines/EVE"
    echo "  # train_VAE + compute_evol_indices, then:"
    echo "  python scripts/prepare_eve_scores.py --csv-path ... --output $EVE_SCORES \\"
    echo "      --data $DATA_DIR --max-seq-len $MAX_SEQ_LEN --ref-from-group 1"
    exit 1
fi

if [[ -z "$OUT_JSON" ]]; then
    OUT_JSON="checkpoints/eval_eve_$(basename "$(dirname "$CHECKPOINT")").json"
fi

echo "Start: $(date)"
echo "checkpoint=$CHECKPOINT"
echo "data=$DATA_DIR  eve=$EVE_SCORES  L=$MAX_SEQ_LEN"
echo "out=$OUT_JSON"

$PYTHON scripts/eval_eve_baseline.py \
    --checkpoint "$CHECKPOINT" \
    --data "$DATA_DIR" \
    --max-seq-len "$MAX_SEQ_LEN" \
    --eve-scores "$EVE_SCORES" \
    --mutation-rate-scale 0.3 \
    --n-steps 100 \
    --max-trees 20 \
    --out "$OUT_JSON"

echo "Done: $(date)  -> $OUT_JSON"
