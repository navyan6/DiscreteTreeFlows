#!/bin/bash
#SBATCH --job-name=covid_train
#SBATCH --partition=b200-mig90
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=16:00:00
#SBATCH --output=logs/covid_train_%j.log
#SBATCH --error=logs/covid_train_%j.log
#
# COVID Spike precompute + bridge-matched training (GPU). Run AFTER
# slurm_covid_pipeline.sh completes (needs group_NNN_rooted.nwk / _anc_aa.fasta
# / _bl.json in each split dir).
#
# Separate checkpoint dir from H3N2 (checkpoints/h3n2_v2) -- does not touch or
# overwrite those weights; trains from scratch (train.py has no warm-start
# mechanism, and Spike's longer max-seq-len changes shape-dependent layers
# anyway). Spike (~1273 aa) is much longer than HA (~566 aa), so --max-seq-len
# must be raised everywhere it's threaded through; ESM2 uses rotary position
# embeddings so it handles the longer length fine (verified with a real
# forward pass locally, no truncation/crash at 1273 aa).
#
# --per-site-pos-emb intentionally OMITTED: on H3N2 it regressed recovery to
# 0 by overfitting to training-era mutation hotspots that didn't generalize
# to held-out data -- no reason to expect it'd do better here.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints/covid_v1

echo "Start: $(date)"

# 1. ESM-2 embeddings + reference (ESM-only) mutation rates, per split
for split in train val test; do
    echo "=== precompute $split ==="
    $PYTHON scripts/precompute_plm.py       --data data/covid/$split
    $PYTHON scripts/precompute_ref_rates.py --data data/covid/$split --max-seq-len 1280
done

# 2. Bridge-matched training on the geographic split
$PYTHON -u scripts/train.py \
    --data        data/covid/train \
    --val-data    data/covid/val \
    --test-data   data/covid/test \
    --max-seq-len 1280 \
    --epochs      100 \
    --patience    50 \
    --bridge-c    1.0 \
    --lambda-mut  5.0 \
    --ckpt-dir    checkpoints/covid_v1

echo "Done: $(date)"
