#!/bin/bash
#SBATCH --job-name=track_a
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/track_a_%j.log
#SBATCH --error=logs/track_a_%j.log
#
# Track A1 Coverage@K (benchmarks/track_a.py) for paper Table 3.
# Default: covid_v5_mutrec @ K ∈ {100,500,1000}. Override via env.
#
#   CKPT_NAME=covid_v5_mutrec KS="100 500 1000" \
#     sbatch --qos=mig-max --job-name=track_a_v5 scripts/slurm_track_a.sh
#   CKPT_NAME=covid_v6_hotspot KS="100 500 1000" \
#     sbatch --qos=mig-max --job-name=track_a_v6 scripts/slurm_track_a.sh
#   # smoke: KS=10 N_GROUPS=2 NSTEPS=20 MAX_LEAVES=50

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-$HF_HOME}"
export TORCH_HOME="${TORCH_HOME:-$LABHOME/torch_cache}"

cd ~/DiscreteTreeFlows
mkdir -p logs benchmarks/results

CKPT_NAME="${CKPT_NAME:-covid_v5_mutrec}"
DATA="${DATA:-data/covid/test}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-1280}"
KS="${KS:-100 500 1000}"
N_GROUPS="${N_GROUPS:-20}"
NSTEPS="${NSTEPS:-50}"
MAX_LEAVES="${MAX_LEAVES:-400}"
BRANCH_RATE_SCALE="${BRANCH_RATE_SCALE:-6.0}"
MUTATION_RATE_SCALE="${MUTATION_RATE_SCALE:-0.04}"
SKIP_EXISTING="${SKIP_EXISTING:-1}"

CKPT="checkpoints/${CKPT_NAME}/best.pt"
if [[ ! -f "$CKPT" && -f "$LABHOME/checkpoints_backup/${CKPT_NAME}/best.pt" ]]; then
  mkdir -p "checkpoints/${CKPT_NAME}"
  cp -n "$LABHOME/checkpoints_backup/${CKPT_NAME}/best.pt" "$CKPT"
fi
if [[ ! -f "$CKPT" ]]; then
  echo "ERROR: missing $CKPT"
  exit 1
fi

echo "Start: $(date)"
echo "ckpt=$CKPT data=$DATA max_seq_len=$MAX_SEQ_LEN KS=[$KS] n_groups=$N_GROUPS"
ls -lah "$CKPT"

for K in $KS; do
  OUT="benchmarks/results/track_a_${CKPT_NAME}_K${K}.json"
  if [[ "$SKIP_EXISTING" == "1" && -f "$OUT" ]]; then
    echo "SKIP K=$K (exists: $OUT)"
    continue
  fi
  echo ""
  echo "############################################################"
  echo "# track_a  $CKPT_NAME  K=$K"
  echo "############################################################"
  $PYTHON benchmarks/track_a.py \
    --checkpoint "$CKPT" \
    --data "$DATA" \
    --K "$K" \
    --n-groups "$N_GROUPS" \
    --n-steps "$NSTEPS" \
    --max-leaves "$MAX_LEAVES" \
    --branch-rate-scale "$BRANCH_RATE_SCALE" \
    --mutation-rate-scale "$MUTATION_RATE_SCALE" \
    --max-seq-len "$MAX_SEQ_LEN" \
    --out "$OUT"
done

echo "Done: $(date)"
echo "Artifacts: benchmarks/results/track_a_${CKPT_NAME}_K*.json"
