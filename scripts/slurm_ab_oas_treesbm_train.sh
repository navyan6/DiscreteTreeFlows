#!/bin/bash
#SBATCH --job-name=ab_oas_train
#SBATCH --partition=b200-mig90
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/ab_oas_train_%j.log
#SBATCH --error=logs/ab_oas_train_%j.log
#
# Train TreeSBM on OAS / Homo_sapiens clone trees (ab_t5_1m → ab_clones_1m).
# Entropy ON + CDR IMGT hotspot mask (NOT pathogen PMC/flu lit).
#
# Distinct from antibody_benchmark ab_dasm_v1 (naive-stripped DASM PCPs).
#
# Usage (after preprocess):
#   sbatch --dependency=afterok:<prep_job> scripts/slurm_ab_oas_treesbm_train.sh
#
# Env:
#   DATA=data/ab_clones_1m CKPT_DIR=checkpoints/ab_oas_1m_v1
#   MUT_HOTSPOT_MASK=results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt

set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
export HF_HOME="${HF_HOME:-$LABHOME/hf_cache}"

cd ~/DiscreteTreeFlows
mkdir -p logs

DATA="${DATA:-data/ab_clones_1m}"
CKPT_DIR="${CKPT_DIR:-checkpoints/ab_oas_1m_v1}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-160}"
LAMBDA_MUT="${LAMBDA_MUT:-12.0}"
LAMBDA_CONS="${LAMBDA_CONS:-0.5}"
LAMBDA_SEMI="${LAMBDA_SEMI:-0.05}"
# Ab genetic BLs ≠ calendar time → disable branch-length / clock loss
LAMBDA_BR="${LAMBDA_BR:-0}"
MUT_NORMALIZE="${MUT_NORMALIZE:-count}"
ENTROPY_ALPHA="${ENTROPY_ALPHA:-3.0}"
ENTROPY_ALPHA_CONS="${ENTROPY_ALPHA_CONS:-1.0}"
MUT_HOTSPOT_MASK="${MUT_HOTSPOT_MASK:-results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt}"
MUT_HOTSPOT_WEIGHT="${MUT_HOTSPOT_WEIGHT:-5.0}"
MUT_HOTSPOT_FORCE="${MUT_HOTSPOT_FORCE:-1}"
PSSM_GATE="${PSSM_GATE:-1}"
MUT_AA_EMB="${MUT_AA_EMB:-1}"
MUT_AA_EMB_DIM="${MUT_AA_EMB_DIM:-16}"
SKIP_PRECOMPUTE="${SKIP_PRECOMPUTE:-0}"
# Gap-filled plm/ref_rates caches from tip-only anc_aa are poisoned — refresh by default.
REFRESH_CACHES="${REFRESH_CACHES:-1}"
RESUME="${RESUME:-0}"

mkdir -p "$CKPT_DIR"

if [[ ! -d "$DATA/train" ]]; then
  echo "ERROR: missing $DATA/train — run slurm_ab_t5_1m_preprocess.sh first" >&2
  exit 1
fi

# Build stub CDR mask if preprocess mask missing
if [[ ! -f "$MUT_HOTSPOT_MASK" ]]; then
  echo "WARN: missing $MUT_HOTSPOT_MASK — building stub"
  $PYTHON scripts/build_ab_cdr_hotspot_mask.py --stub --max-seq-len "$MAX_SEQ_LEN" \
    --out "$MUT_HOTSPOT_MASK"
fi

EXTRA_ARGS=()
if [[ "$PSSM_GATE" == "1" || "$PSSM_GATE" == "true" ]]; then
  EXTRA_ARGS+=(--pssm-gate)
fi
if [[ "$MUT_AA_EMB" == "1" || "$MUT_AA_EMB" == "true" ]]; then
  EXTRA_ARGS+=(--mut-aa-emb --mut-aa-emb-dim "$MUT_AA_EMB_DIM")
fi
EXTRA_ARGS+=(--entropy-weight-alpha-cons "$ENTROPY_ALPHA_CONS")
EXTRA_ARGS+=(--mut-hotspot-mask "$MUT_HOTSPOT_MASK")
EXTRA_ARGS+=(--mut-hotspot-weight "$MUT_HOTSPOT_WEIGHT")
if [[ "$MUT_HOTSPOT_FORCE" == "1" || "$MUT_HOTSPOT_FORCE" == "true" ]]; then
  EXTRA_ARGS+=(--mut-hotspot-force)
fi
if [[ "$RESUME" == "1" || "$RESUME" == "true" ]]; then
  EXTRA_ARGS+=(--resume)
fi

echo "=== Ab OAS TreeSBM train ==="
echo "host=$(hostname) date=$(date -Is) job=${SLURM_JOB_ID:-local}"
echo "data=$DATA ckpt=$CKPT_DIR max_seq_len=$MAX_SEQ_LEN"
echo "entropy ON (empirical) + CDR mask=$MUT_HOTSPOT_MASK"
echo "lambda_mut=$LAMBDA_MUT lambda_cons=$LAMBDA_CONS lambda_semi=$LAMBDA_SEMI lambda_br=$LAMBDA_BR"
echo "resume=$RESUME"
echo "REFRESH_CACHES=$REFRESH_CACHES SKIP_PRECOMPUTE=$SKIP_PRECOMPUTE"
ls -la "$MUT_HOTSPOT_MASK"
"$PYTHON" -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"

# Sanity: TreeDataset must see complete groups before precompute/train
# (heredoc — do not interpolate paths inside python -c quotes)
for split in train val; do
  n=$(DATA_SPLIT="$DATA/$split" MAX_SEQ_LEN="$MAX_SEQ_LEN" $PYTHON - <<'PY'
import os, sys
from src.dataset import TreeDataset
ds = TreeDataset(os.environ["DATA_SPLIT"], max_seq_len=int(os.environ["MAX_SEQ_LEN"]))
item = ds[0]
print(f"smoke_group {item['group']} n_nodes={len(item['node_ids'])}", file=sys.stderr, flush=True)
# no all-gap internals after parent fill
n_gap = sum(1 for s in item["seqs"].values() if s and set(s) <= set("-."))
print(f"smoke_allgap_nodes={n_gap}", file=sys.stderr, flush=True)
if n_gap:
    raise SystemExit(f"ERROR: {n_gap} all-gap node seqs after fill_missing_node_seqs")
print(len(ds.groups))
PY
)
  n=$(printf '%s\n' "$n" | tail -1)
  echo "TreeDataset $split: Found $n complete groups"
  if [[ "$n" -le 0 ]]; then
    echo "ERROR: TreeDataset found 0 trees in $DATA/$split — run ab_postprocess_clones_for_treesbm.py first" >&2
    exit 1
  fi
done

if [[ "$REFRESH_CACHES" == "1" || "$REFRESH_CACHES" == "true" ]]; then
  echo "=== refreshing plm/ref_rates caches (parent-fill seqs) ==="
  for split in train val; do
    find "$DATA/$split" -maxdepth 1 -name 'group_*_plm.pt' -delete
    find "$DATA/$split" -maxdepth 1 -name 'group_*_ref_rates*.pt' -delete
  done
fi

if [[ "$SKIP_PRECOMPUTE" != "1" ]]; then
  for split in train val; do
    echo "=== precompute_plm $split ==="
    $PYTHON -u scripts/precompute_plm.py --data "$DATA/$split" --overwrite
    echo "=== precompute_ref_rates $split ==="
    $PYTHON -u scripts/precompute_ref_rates.py \
      --data "$DATA/$split" \
      --max-seq-len "$MAX_SEQ_LEN" \
      --overwrite
  done
else
  echo "SKIP_PRECOMPUTE=1 — using existing *_plm.pt / *_ref_rates.pt"
fi

echo "=== train.py ==="
$PYTHON -u scripts/train.py \
  --data "$DATA/train" \
  --val-data "$DATA/val" \
  --test-data "$DATA/test" \
  --max-seq-len "$MAX_SEQ_LEN" \
  --epochs 150 \
  --patience 40 \
  --lr 1e-4 \
  --bridge-c 1.0 \
  --lambda-mut "$LAMBDA_MUT" \
  --lambda-cons "$LAMBDA_CONS" \
  --lambda-top 0.5 \
  --lambda-br "$LAMBDA_BR" \
  --lambda-semi "$LAMBDA_SEMI" \
  --mut-normalize "$MUT_NORMALIZE" \
  --per-site-pos-emb \
  --use-site-entropy \
  --use-entropy-loss-weighting \
  --use-entropy-cons-weighting \
  --entropy-source empirical \
  --entropy-weight-alpha "$ENTROPY_ALPHA" \
  --entropy-weight-floor 1.0 \
  --n-t-samples 4 \
  --ckpt-dir "$CKPT_DIR" \
  "${EXTRA_ARGS[@]}"

echo "=== done ==="
ls -la "$CKPT_DIR/"
echo "Done: $(date -Is)"
echo "Eval plan:"
echo "  1) Leaf/maturation on $DATA/test (scripts/eval_ab_maturation.py / slurm_ab_t5_eval.sh)"
echo "  2) Keep Rodriguez 82 antibody_benchmark as separate CoSiNE/DASM table — do not replace"
echo "  3) Optional: rollout antibody_benchmark only if ckpt shows non-collapsed mut/rate/cons"
