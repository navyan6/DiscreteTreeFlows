#!/bin/bash
# Pull SPLIT_PROTOCOL.json + group CSVs from Betty (no FASTA). Run from laptop.
set -euo pipefail
REMOTE="${BETTY:-nnori@login.betty.parcc.upenn.edu}"
RROOT="${RROOT:-~/DiscreteTreeFlows}"
LOCAL="${LOCAL:-$(cd "$(dirname "$0")/.." && pwd)}"

DIRS=(
  data/h3n2
  data/h1n1
  data/covid
  data/hiv_geo
  data/hiv_temporal
  data/covid_epidemic
  data/h3n2_epidemic
  data/h1n1_epidemic
  data/filo_l
)

for d in "${DIRS[@]}"; do
  echo "=== $d ==="
  ssh "$REMOTE" "test -f $RROOT/$d/SPLIT_PROTOCOL.json" || { echo "  MISSING protocol"; continue; }
  mkdir -p "$LOCAL/$d/train" "$LOCAL/$d/val" "$LOCAL/$d/test"
  rsync -avz \
    "$REMOTE:$RROOT/$d/SPLIT_PROTOCOL.json" \
    "$LOCAL/$d/"
  for split in train val test; do
    rsync -avz --include='*_group_*.csv' --exclude='*' \
      "$REMOTE:$RROOT/$d/$split/" "$LOCAL/$d/$split/" 2>/dev/null || true
  done
done

rsync -avz "$REMOTE:$RROOT/benchmarks/results/tables/epidemic_split_audit.json" \
  "$LOCAL/benchmarks/results/tables/" 2>/dev/null || true

echo "Done. Review with: git status data/"
