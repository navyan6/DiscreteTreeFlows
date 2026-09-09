#!/bin/bash
#SBATCH --job-name=eval_covid
#SBATCH --partition=b200-mig45
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/eval_covid_%j.log
#SBATCH --error=logs/eval_covid_%j.log
#
# Covid-only EVEscape enrichment eval (GPU): covid_v1 (baseline) vs
# covid_v2_entropy (empirical entropy), both WITH EVEscape RBD scoring.
# h1n1 already ran successfully, so not repeated here. --mutation-rate-scale
# 0.3 per the smoke-test finding.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs
echo "Start: $(date)"

for ckpt in covid_v1 covid_v2_entropy; do
    echo ""
    echo "############################################################"
    echo "# $ckpt"
    echo "############################################################"
    $PYTHON scripts/eval_evescape_enrichment.py \
        --checkpoint checkpoints/$ckpt/best.pt \
        --data data/covid/test \
        --max-seq-len 1280 \
        --evescape data/covid/evescape_spike_rbd.pt \
        --mutation-rate-scale 0.3 \
        --n-steps 100 \
        --max-trees 20
done

echo ""
echo "Done: $(date)"
