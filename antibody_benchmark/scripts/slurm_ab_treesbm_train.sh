#!/bin/bash
#SBATCH --job-name=ab_treesbm_train
#SBATCH --partition=b200-mig90
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=/vast/projects/pranam/lab/nnori/antibody_benchmark/logs/ab_treesbm_train_%j.out
#SBATCH --error=/vast/projects/pranam/lab/nnori/antibody_benchmark/logs/ab_treesbm_train_%j.err
#
# Train TreeSBM on DASM heavy Ab trees (Tang+VanWinkle; NOT Rodriguez held-out).
# After success: checkpoints/ab_dasm_v1/best.pt

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:${PATH:-}
export AB_DASM_DIR="${AB_DASM_DIR:-$HOME/antibody_benchmark_raw/dasm/extracted/dasm-experiments-data}"
LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
PY="${PY:-/vast/projects/pranam/lab/nnori/.conda/envs/treesbm/bin/python}"

cd "$HOME/DiscreteTreeFlows"
mkdir -p "$LABHOME/antibody_benchmark/logs" logs checkpoints/ab_dasm_v1 data/ab_dasm_trees

echo "=== Ab TreeSBM train ==="
echo "host=$(hostname) date=$(date -Is) job=$SLURM_JOB_ID"
echo "GPU=$CUDA_VISIBLE_DEVICES AB_DASM_DIR=$AB_DASM_DIR"
"$PY" -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"

# 1) Export DASM → TreeDataset (skip if already present with meta)
if [[ ! -f data/ab_dasm_trees/export_meta.json ]]; then
  echo "=== export train trees ==="
  "$PY" -u antibody_benchmark/scripts/export_treesbm_train_data.py \
    --dasm-dir "$AB_DASM_DIR" \
    --out data/ab_dasm_trees \
    --max-trees 800 \
    --val-frac 0.1 \
    --seed 42
else
  echo "=== reuse existing data/ab_dasm_trees ==="
  cat data/ab_dasm_trees/export_meta.json
fi

MAX_SEQ_LEN="${MAX_SEQ_LEN:-200}"
echo "max_seq_len=$MAX_SEQ_LEN"

# 2) Precompute ESM embeddings + ref rates
for split in train val; do
  echo "=== precompute_plm $split ==="
  "$PY" -u scripts/precompute_plm.py --data "data/ab_dasm_trees/$split"
  echo "=== precompute_ref_rates $split ==="
  "$PY" -u scripts/precompute_ref_rates.py \
    --data "data/ab_dasm_trees/$split" \
    --max-seq-len "$MAX_SEQ_LEN"
done

# 3) Bridge-matched training (hparams mirror scripts/slurm_h3n2_train.sh)
echo "=== train.py ==="
"$PY" -u scripts/train.py \
  --data data/ab_dasm_trees/train \
  --val-data data/ab_dasm_trees/val \
  --epochs 150 \
  --patience 40 \
  --lr 1e-4 \
  --bridge-c 1.0 \
  --lambda-mut 5.0 \
  --lambda-top 0.5 \
  --lambda-br 0.1 \
  --per-site-pos-emb \
  --max-seq-len "$MAX_SEQ_LEN" \
  --n-t-samples 4 \
  --ckpt-dir checkpoints/ab_dasm_v1

echo "=== done ==="
ls -la checkpoints/ab_dasm_v1/
echo "Done: $(date -Is)"
