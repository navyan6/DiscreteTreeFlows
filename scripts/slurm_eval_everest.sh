#!/bin/bash
# EVEREST forecasting eval on epidemic + Wave-1 checkpoints (Betty).
set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs benchmarks/results/everest external

if [[ ! -d external/priority-viruses/data/forecasting_dataset ]]; then
  echo "Cloning priority-viruses forecasting_dataset..."
  git clone --depth 1 https://github.com/debbiemarkslab/priority-viruses external/priority-viruses
fi

run_eval() {
  local name="$1" ckpt="$2" data="$3" virus="$4" maxlen="$5" everest_sub="$6"
  local out="benchmarks/results/everest/${name}.json"
  python scripts/eval_everest_forecasting.py \
    --checkpoint "$ckpt" \
    --data "$data" \
    --virus "$virus" \
    --max-seq-len "$maxlen" \
    --everest-root "external/priority-viruses/data/forecasting_dataset/${everest_sub}" \
    --K 16 --max-trees 20 \
    --out "$out"
}

# Epidemic ckpts (after retrain completes)
run_eval h3n2_v4_epidemic \
  checkpoints/h3n2_v4_epidemic_mutrec/best.pt \
  data/h3n2_epidemic/test h3n2 566 H3N2 || true

run_eval h1n1_v3_epidemic \
  checkpoints/h1n1_v3_epidemic_mutrec/best.pt \
  data/h1n1_epidemic/test h1n1 566 H1N1 || true

run_eval covid_v6_epidemic \
  checkpoints/covid_v6_epidemic_mutrec/best.pt \
  data/covid_epidemic/test covid 1280 SARS-CoV-2 || true

# Wave-1 ablation (geo/temporal paper splits)
run_eval h3n2_v3_wave1 \
  checkpoints/h3n2_v3_lit_hotspot/best.pt \
  data/h3n2/test h3n2 566 H3N2 || true

run_eval h1n1_wave1 \
  checkpoints/h1n1_v2_mutrec/best.pt \
  data/h1n1/test h1n1 566 H1N1 || true

run_eval covid_v5_wave1 \
  checkpoints/covid_v5_mutrec/best.pt \
  data/covid/test covid 1280 SARS-CoV-2 || true

echo "EVEREST eval batch done — see benchmarks/results/everest/"
