#!/bin/bash
#SBATCH --job-name=phylovae
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/vast/home/n/nnori/DiscreteTreeFlows/logs/phylovae_%j.log
#SBATCH --error=/vast/home/n/nnori/DiscreteTreeFlows/logs/phylovae_%j.log
#
# PhyloVAE: train a topology-prior model for ONE N on our anonymized
# train-set topology pool, then sample a fresh pool from it. CPU-only by the
# repo's own design (its README calls it "a light CPU-based implementation").
# Own conda env (never treesbm -- see benchmarks/EXTERNAL.md), lives outside
# DiscreteTreeFlows.
#
# One N per job, not all three sequentially: N=16/32/64 have very different
# per-iteration cost (bigger trees = bigger tensors -- observed 15.22/9.81/
# 11.25 it/s respectively), so packing all three into one job risks the
# smaller N's eating the shared time budget before the largest N gets enough
# of it. This happened: N=64 only reached 45k/200k iters (22%) in a 12h job
# that N=16+N=32 alone had already used 9h19m of. Submit once per N instead;
# they can run concurrently and each gets its own full budget.
#
# Usage: sbatch scripts/slurm_phylovae.sh 16   (or 32, or 64)
#
# Prereq:
#   python benchmarks/heldout/export_train_topologies.py --N 16 32 64
#   cp benchmarks/external_adapters/phylovae_sample.py \
#       $REPO_DIR/phylovae_sample.py
#   conda env create -f $REPO_DIR/environment.yaml   (one-time)

set -e
N="${1:?usage: sbatch slurm_phylovae.sh <N>  (16, 32, or 64)}"
REPO_DIR=/vast/projects/pranam/lab/nnori/baselines/PhyloVAE
POOL_SRC=/vast/home/n/nnori/DiscreteTreeFlows/benchmarks/external_pools
POOL_OUT=/vast/home/n/nnori/DiscreteTreeFlows/benchmarks/external_pools/sampled
mkdir -p "$POOL_OUT" /vast/home/n/nnori/DiscreteTreeFlows/logs

source activate phylovae
cd "$REPO_DIR"

DATASET="treesbm_N${N}"
echo "=== N=${N}: placing .trprobs ==="
mkdir -p "data/short_run_data/${DATASET}/rep_1"
cp "${POOL_SRC}/train_topologies_N${N}.trprobs" \
   "data/short_run_data/${DATASET}/rep_1/${DATASET}.trprobs"

echo "=== N=${N}: process_data ==="
python -c "from src.datasets import process_data; process_data('${DATASET}', 1)"

echo "=== N=${N}: train ==="
python main.py base.mode=train data.dataset="${DATASET}" data.rep_id=1 \
    decoder.num_layers=4 decoder.latent_dim=2 \
    objective.batch_size=10 objective.n_particles=32

echo "=== N=${N}: sample ==="
CKPT=$(find "${REPO_DIR}/results/tde/${DATASET}/rep_1" -name final.pt | head -1)
if [ -z "$CKPT" ]; then
    echo "ERROR: no final.pt found for N=${N} under results/tde/${DATASET}/rep_1"
    exit 1
fi
python phylovae_sample.py \
    base.mode=train data.dataset="${DATASET}" data.rep_id=1 \
    decoder.num_layers=4 decoder.latent_dim=2 \
    objective.batch_size=10 objective.n_particles=32 \
    --checkpoint "${CKPT}" --ntips "${N}" \
    --n-samples 300 --out "${POOL_OUT}/phylovae_N${N}.nwk"

echo "Done: $(date)"
