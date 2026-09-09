#!/bin/bash
# Local Track A prep: screen Brazil + origin groups, assemble cases, smoke-eval Gamma.
# Safe to run without GPU. Betty: then bash scripts/betty_submit_voc_threat_panel.sh

set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/voc_threat_panel

echo "=== 1. Brazil geo-test screen ==="
python scripts/screen_voc_threat_trees.py \
  --data-dir data/covid/test --split test \
  --vocs Gamma,Zeta,Alpha,Beta,Delta,Omicron_BA1,Lambda,Mu \
  --out results/voc_threat_panel/brazil_screen.json \
  --candidates-out results/voc_threat_panel/candidates_test.json

echo "=== 2. Origin-country train screen ==="
python scripts/screen_voc_threat_trees.py \
  --data-dir data/covid/train --split train \
  --groups 34,289,188,189,190,191,192,193,300,301,302,303,304,308,309,323 \
  --vocs Gamma,Zeta,Alpha,Beta,Delta,Omicron_BA1,Lambda,Mu,Kappa \
  --out results/voc_threat_panel/origin_screen.json \
  --candidates-out results/voc_threat_panel/candidates_train.json

echo "=== 3. Assemble cases ==="
python scripts/prepare_voc_origin_cases.py \
  --candidates results/voc_threat_panel/candidates_test.json \
  --candidates results/voc_threat_panel/candidates_train.json \
  --case Gamma:test:3 \
  --case Beta:train:302 \
  --case Delta:train:191 \
  --case Omicron_BA1:train:309 \
  --case Alpha:train:323 \
  --case Lambda:train:289 \
  --case Mu:train:34 \
  --max-per-voc 2

echo "=== 4. Smoke-eval existing Gamma gens (if present) ==="
GEN=results/covid_tree_viz/group_003_generated_matched.fasta
OBS=results/covid_tree_viz/group_003_observed_anc_aa.fasta
if [[ -f "$GEN" && -f "$OBS" ]]; then
  python scripts/eval_voc_threat_recovery.py \
    --voc Gamma \
    --gen-fasta "$GEN" \
    --obs-fasta "$OBS" \
    --ref-tip MZ397166.1 \
    --topk 10 \
    --out results/voc_threat_panel/eval_gamma_g003_existing.json
else
  echo "SKIP smoke-eval (missing $GEN or $OBS)"
fi

echo "=== Done. On Betty: DRY=1 bash scripts/betty_submit_voc_threat_panel.sh ==="
ls -la results/voc_threat_panel/cases_manifest.json results/voc_threat_panel/candidates_*.json
