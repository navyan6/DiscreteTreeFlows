#!/bin/bash
#SBATCH --job-name=baselines
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/baselines_%j.log
#SBATCH --error=logs/baselines_%j.log
#
# Native baseline table (Table 2 style): NeutralBD / EmpiricalBD / PLMPrior,
# plus TreeSBM when a checkpoint is passed. Adapted ARTreeFormer / PhyloVAE and
# PhylaFlow native (`phylaflow`) / adapted rows are added by run_table.py when
# benchmarks/external_pools/sampled/*.nwk exist (see EXTERNAL.md).
#
# Walltime: 24h (Wave-1 job 7323289 timed out at 8h before TreeSBM landed in
# table_empirical_main.csv). Default QOS mig-max avoids MaxGRESPerAccount on qos=mig.
#
# Quartet column needs tqdist (requirements.txt / treesbm env). Missing -> NaN,
# never faked. This preamble verifies/installs tqdist if absent.
#
# Usage / resubmit (always pass --qos=mig-max if overriding header):
#   # BD + PLM only (no TreeSBM row):
#   sbatch --qos=mig-max scripts/slurm_baselines.sh
#   # Include TreeSBM (H3N2) — primary Table 2 fill:
#   sbatch --qos=mig-max scripts/slurm_baselines.sh checkpoints/h3n2_v2/best.pt
#   # Explicit time override (same as header):
#   sbatch --qos=mig-max --time=24:00:00 scripts/slurm_baselines.sh checkpoints/h3n2_v2/best.pt
#   # COVID track:
#   DATA=data/covid/test TRAIN_DATA=data/covid/train \
#     OUT=benchmarks/results/results_baselines_covid.csv \
#     sbatch --qos=mig-max scripts/slurm_baselines.sh checkpoints/covid_v3_cons/best.pt
#
# Env overrides: CKPT, DATA/TEST_DATA, TRAIN_DATA, OUT, PARAMS, N_LIST, K, M, MAX_ROOTS

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"

cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results

# Positional arg wins; else CKPT env; else empty (BD+PLM only).
CKPT="${1:-${CKPT:-}}"
TEST_DATA="${DATA:-${TEST_DATA:-data/h3n2/test}}"
TRAIN_DATA="${TRAIN_DATA:-data/h3n2/train}"
PARAMS="${PARAMS:-benchmarks/results/params.json}"
OUT="${OUT:-benchmarks/results/results_baselines.csv}"
TABLE_OUT="${TABLE_OUT:-benchmarks/results/tables}"
N_LIST="${N_LIST:-16 32 64}"
K="${K:-20}"
M="${M:-20}"
MAX_ROOTS="${MAX_ROOTS:-100}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"
# Space-separated regimes for sim track; set REGIMES="" to skip (empirical-only KL).
REGIMES="${REGIMES-neutral JTT WAG LG}"

resolve_ckpt() {
    local c="$1"
    [ -z "$c" ] && return 0
    if [ -f "$c" ]; then echo "$c"; return; fi
    local name
    name=$(basename "$(dirname "$c")")
    local backup="$LABHOME/checkpoints_backup/${name}/best.pt"
    if [ -f "$backup" ]; then
        mkdir -p "$(dirname "$c")"
        cp -n "$backup" "$c" 2>/dev/null || cp "$backup" "$c"
        echo "$c"
        return
    fi
    if [ -f "$backup" ]; then echo "$backup"; return; fi
    echo "$c"
}

echo "Start: $(date)"
echo "Host: $(hostname)  PWD: $(pwd)"
echo "test=$TEST_DATA train=$TRAIN_DATA out=$OUT N=($N_LIST) K=$K M=$M max-roots=$MAX_ROOTS max_seq_len=$MAX_SEQ_LEN regimes=($REGIMES)"

# tqdist: required for Quartet; install if missing (do not fake NaNs away).
if ! $PYTHON -c "import tqdist" 2>/dev/null; then
    echo "tqdist missing — pip installing into treesbm env"
    $PYTHON -m pip install -q 'tqdist>=1.0.0'
fi
$PYTHON -c "import tqdist; print('tqdist OK', tqdist.__file__)"

if [ ! -f "$PARAMS" ]; then
    echo "Fitting BD/subst params on $TRAIN_DATA -> $PARAMS"
    $PYTHON benchmarks/fit_params.py --train-data "$TRAIN_DATA" --out "$PARAMS"
fi

EXTRA=()
if [ -n "$CKPT" ]; then
    CKPT_RESOLVED="$(resolve_ckpt "$CKPT")"
    if [ -f "$CKPT_RESOLVED" ]; then
        echo "TreeSBM checkpoint: $CKPT_RESOLVED"
        EXTRA+=(--checkpoint "$CKPT_RESOLVED")
    else
        echo "ERROR: checkpoint not found: $CKPT (tried lab backup under $LABHOME/checkpoints_backup/)"
        exit 1
    fi
else
    echo "No checkpoint arg/CKPT — NeutralBD + EmpiricalBD + PLMPrior only (TreeSBM skipped)"
fi

# shellcheck disable=SC2086
$PYTHON benchmarks/run_table.py \
    --test-data "$TEST_DATA" \
    --train-data "$TRAIN_DATA" \
    --params "$PARAMS" \
    --N $N_LIST --K "$K" --M "$M" --max-roots "$MAX_ROOTS" \
    --max-seq-len "$MAX_SEQ_LEN" \
    --regimes $REGIMES \
    --out "$OUT" \
    "${EXTRA[@]}"

echo "Aggregating with make_table.py ..."
$PYTHON benchmarks/make_table.py --results "$OUT" --out "$TABLE_OUT" --ci se

echo "Done: $(date)"
echo "Tables under $TABLE_OUT/"
