#!/bin/bash
#SBATCH --job-name=eval_suite
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/eval_suite_%j.log
#SBATCH --error=logs/eval_suite_%j.log
#
# Aggregate eval across the three just-completed runs: covid baseline
# (covid_v1), covid + empirical entropy (covid_v2_entropy, WITH EVEscape
# RBD enrichment), and h1n1 (h1n1_v1, no EVEscape -- matrix is COVID-spike-
# RBD only). --mutation-rate-scale 0.3 per the smoke-test finding (rate=1.0
# over-mutates: retention 0.65 vs 0.88 at 0.3).
#
# Each run writes checkpoints/eval_enrichment_<ckpt-dir-name>.json + prints
# a summary block to this log. max-trees 20 for a same-day turnaround; raise
# later for a fuller number if there's time.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Start: $(date)"

echo ""
echo "############################################################"
echo "# 1/3  covid_v1 (baseline, no entropy)"
echo "############################################################"
$PYTHON scripts/eval_evescape_enrichment.py \
    --checkpoint checkpoints/covid_v1/best.pt \
    --data data/covid/test \
    --max-seq-len 1280 \
    --evescape data/covid/evescape_spike_rbd.pt \
    --mutation-rate-scale 0.3 \
    --n-steps 100 \
    --max-trees 20

echo ""
echo "############################################################"
echo "# 2/3  covid_v2_entropy (empirical column entropy)"
echo "############################################################"
$PYTHON scripts/eval_evescape_enrichment.py \
    --checkpoint checkpoints/covid_v2_entropy/best.pt \
    --data data/covid/test \
    --max-seq-len 1280 \
    --evescape data/covid/evescape_spike_rbd.pt \
    --mutation-rate-scale 0.3 \
    --n-steps 100 \
    --max-trees 20

echo ""
echo "############################################################"
echo "# 3/3  h1n1_v1 (empirical column entropy, no EVEscape -- RBD-only matrix)"
echo "############################################################"
$PYTHON scripts/eval_evescape_enrichment.py \
    --checkpoint checkpoints/h1n1_v1/best.pt \
    --data data/h1n1/test \
    --max-seq-len 566 \
    --mutation-rate-scale 0.3 \
    --n-steps 100 \
    --max-trees 20

echo ""
echo "Done: $(date)"
