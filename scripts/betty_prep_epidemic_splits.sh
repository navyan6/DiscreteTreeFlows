#!/bin/bash
# Prepare epidemic v2 split dirs for COVID + flu (no SLURM train). Run on Betty after data exists.
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
cd ~/DiscreteTreeFlows

LAB=/vast/projects/pranam/lab/nnori/DiscreteTreeFlows_data
mkdir -p data/covid data/h1n1
ln -sfn "$LAB/covid/train" data/covid/train
ln -sfn "$LAB/h1n1/train" data/h1n1/train

echo "=== COVID epidemic ==="
if ls data/covid/train/*_spike.fasta &>/dev/null; then
  $PYTHON scripts/prepare_covid_epidemic.py
else
  echo "SKIP COVID: no spike FASTAs"
fi

echo "=== H3N2 epidemic (from lab temporal pools) ==="
if [[ -f "$LAB/h3n2/train/h3n2train.fasta" ]]; then
  $PYTHON scripts/prepare_h3n2_epidemic.py --pool-base "$LAB/h3n2"
else
  $PYTHON scripts/prepare_h3n2_epidemic.py || echo "WARN: h3n2 source FASTAs missing"
fi

echo "=== H1N1 epidemic ==="
if [[ -f "$LAB/h1n1_temporal/train/h1n1ttrain.fasta" ]]; then
  $PYTHON scripts/prepare_h1n1_epidemic.py --pool-base "$LAB/h1n1_temporal"
elif ls data/h1n1/train/*.fasta &>/dev/null; then
  $PYTHON scripts/prepare_h1n1_epidemic.py
else
  echo "WARN: h1n1 source FASTAs missing"
fi

$PYTHON scripts/audit_epidemic_splits.py \
  --data data/covid_epidemic --data data/h3n2_epidemic --data data/h1n1_epidemic \
  --out benchmarks/results/tables/epidemic_split_audit.json

echo "Epidemic prep done: $(date)"
