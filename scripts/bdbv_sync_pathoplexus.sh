#!/bin/bash
# Import Pathoplexus BDBV export as golden test (Bdbv-2026 only; never train).
#
# Usage:
#   bash scripts/bdbv_sync_pathoplexus.sh
#   bash scripts/bdbv_sync_pathoplexus.sh data/bdbv/ebola-bdbv_metadata_*.tsv.gz
#
# Expects paired metadata + aligned-nuc FASTA under data/bdbv/ (or PATHOPLEXUS_SRC).

set -euo pipefail
export PATH="${PATH:-}"
PYTHON="${PYTHON:-python3}"
cd "$(dirname "$0")/.."

SEARCH="${PATHOPLEXUS_SRC:-data/bdbv}"
META="${1:-}"
FASTA="${2:-}"

ARGS=(--search-dir "$SEARCH")
if [[ -n "$META" ]]; then ARGS+=(--metadata "$META"); fi
if [[ -n "$FASTA" ]]; then ARGS+=(--fasta "$FASTA"); fi

$PYTHON scripts/bdbv_ingest_pathoplexus_export.py "${ARGS[@]}"
echo "Golden test ready under data/filo_l/test/ — re-run eval only (no retrain required)."
