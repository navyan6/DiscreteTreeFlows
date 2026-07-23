#!/bin/bash
#SBATCH --job-name=artreeformer
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/vast/home/n/nnori/DiscreteTreeFlows/logs/artreeformer_%j.log
#SBATCH --error=/vast/home/n/nnori/DiscreteTreeFlows/logs/artreeformer_%j.log
#
# ARTreeFormer TDE: train a topology-prior model per N on our anonymized
# train-set topology pool, then sample a fresh pool from it. Own conda env
# (never treesbm -- see benchmarks/EXTERNAL.md), lives outside DiscreteTreeFlows.
#
# Prereq:
#   python benchmarks/heldout/export_train_topologies.py --N 16 32 64
#   cp benchmarks/external_adapters/artreeformer_sample.py \
#       $REPO_DIR/TDE/artreeformer_sample.py
#   conda env create -f $REPO_DIR/environment.yaml   (one-time)

set -e
REPO_DIR=/vast/projects/pranam/lab/nnori/baselines/ARTreeFormer
POOL_SRC=/vast/home/n/nnori/DiscreteTreeFlows/benchmarks/external_pools
POOL_OUT=/vast/home/n/nnori/DiscreteTreeFlows/benchmarks/external_pools/sampled
mkdir -p "$POOL_OUT" /vast/home/n/nnori/DiscreteTreeFlows/logs

source activate artreeformer
cd "$REPO_DIR"

for N in 16 32 64; do
    DATASET="treesbm_N${N}"
    echo "=== N=${N}: placing .trprobs ==="
    mkdir -p "data/short_run_data_DS1-8/${DATASET}/rep_1"
    cp "${POOL_SRC}/train_topologies_N${N}.trprobs" \
       "data/short_run_data_DS1-8/${DATASET}/rep_1/${DATASET}.trprobs"

    echo "=== N=${N}: process_data ==="
    python -c "from datasets import process_data; process_data('${DATASET}', 1)"

    echo "=== N=${N}: train ==="
    cd TDE
    python main.py data.dataset="${DATASET}" data.repo=1 base.mode=train
    cd ..

    echo "=== N=${N}: sample ==="
    CKPT=$(find "${REPO_DIR}/TDE/results/${DATASET}/repo1" -name final.pt | head -1)
    if [ -z "$CKPT" ]; then
        echo "ERROR: no final.pt found for N=${N} under TDE/results/${DATASET}/repo1"
        exit 1
    fi
    cd TDE
    python artreeformer_sample.py --checkpoint "${CKPT}" --ntips "${N}" \
        --n-samples 300 --out "${POOL_OUT}/artreeformer_N${N}.nwk"
    cd ..
done

echo "Done: $(date)"
