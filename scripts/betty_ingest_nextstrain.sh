#!/bin/bash
# Re-run Nextstrain Charon metadata merge + L extract + outbreak splits (Betty login node).
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
cd ~/DiscreteTreeFlows
mkdir -p logs data/bdbv/references results/bdbv_l_mask

GROUP_SIZE="${GROUP_SIZE:-80}"
MIN_GROUP="${MIN_GROUP:-5}"

echo "=== Nextstrain Charon metadata ==="
$PYTHON scripts/download_ebolavirus_nextstrain.py

echo "=== L extract + window + lit mask (GenBank efetch for new accessions) ==="
$PYTHON scripts/bdbv_extract_l.py
$PYTHON scripts/bdbv_define_l_window.py
$PYTHON scripts/build_bdbv_l_lit_hotspot_mask.py

echo "=== Outbreak splits (group_size=$GROUP_SIZE min_group=$MIN_GROUP) ==="
$PYTHON scripts/prepare_filo_outbreak.py \
  --out-base data/filo_l --group-size "$GROUP_SIZE" --min-group "$MIN_GROUP"
$PYTHON scripts/prepare_filo_outbreak.py --track-b --test-outbreak ebov_wa_2013_2016 \
  --out-base data/filo_l_track_b --group-size "$GROUP_SIZE" --min-group "$MIN_GROUP" || true
$PYTHON scripts/audit_filo_splits.py --data data/filo_l

echo "Nextstrain filo re-ingest done: $(date)"
