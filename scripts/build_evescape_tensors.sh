#!/usr/bin/env bash
# Rebuild EVEscape lookup tensors from Marks official summary CSVs.
# No PDB required. Keeps official log-scale (does NOT z-score).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
PY="${PYTHON:-python3}"
OFF="$ROOT/data/evescape_official"
mkdir -p "$OFF"

echo "Downloading Marks summaries_with_scores (if needed)..."
[[ -f "$OFF/flu_h1_evescape.csv" ]] || curl -fsSL -o "$OFF/flu_h1_evescape.csv" \
  https://raw.githubusercontent.com/OATML-Markslab/EVEscape/main/results/summaries_with_scores/flu_h1_evescape.csv
[[ -f "$OFF/spike_rbd_evescape.csv" ]] || curl -fsSL -o "$OFF/spike_rbd_evescape.csv" \
  https://raw.githubusercontent.com/OATML-Markslab/EVEscape/main/results/summaries_with_scores/spike_rbd_evescape.csv

if [[ -f data/evescape_h1n1_ref_seq.txt ]]; then
  echo "Building data/evescape_h1n1_ha.pt ..."
  $PY scripts/prepare_evescape.py \
    --csv-path "$OFF/flu_h1_evescape.csv" \
    --output data/evescape_h1n1_ha.pt \
    --max-seq-len 566 --pathogen flu_h1 \
    --ref-seq-file data/evescape_h1n1_ref_seq.txt \
    --no-standardize --min-match-rate 0.75
else
  echo "SKIP H1: missing data/evescape_h1n1_ref_seq.txt (or pass --data/--ref-seq manually)"
fi

if [[ -f data/covid/evescape_ref_seq.txt ]]; then
  echo "Building data/covid/evescape_spike_rbd.pt ..."
  $PY scripts/prepare_evescape.py \
    --csv-path "$OFF/spike_rbd_evescape.csv" \
    --output data/covid/evescape_spike_rbd.pt \
    --max-seq-len 1280 --pathogen covid_spike_rbd \
    --ref-seq-file data/covid/evescape_ref_seq.txt \
    --no-standardize
else
  echo "SKIP COVID: missing data/covid/evescape_ref_seq.txt"
fi

echo "Done. Official EVEscape mean should be ≈ −2.3 (not ≈0 z-score, not ≈−8 EVE)."
