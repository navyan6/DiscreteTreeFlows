#!/bin/bash
# Pull formed training trees (nwk + anc_aa + bl.json) from Betty for NCBI-safe splits.
# Usage:
#   bash scripts/pull_formed_trees_from_betty.sh h1n1
#   bash scripts/pull_formed_trees_from_betty.sh h1n1 hiv_geo flub
set -euo pipefail
HOST="${BETTY_HOST:-nnori@login.betty.parcc.upenn.edu}"
REMOTE="${BETTY_ROOT:-~/DiscreteTreeFlows}"
LOCAL="${ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"

pull_split() {
  local name="$1"
  for split in train val test; do
    local rdir="$REMOTE/data/$name/$split"
    local ldir="$LOCAL/data/$name/$split"
    mkdir -p "$ldir"
    echo "=== $name/$split ==="
    rsync -avz --include='group_*_rooted.nwk' --include='group_*_anc_aa.fasta' \
      --include='group_*_bl.json' --include='SPLIT_PROTOCOL.json' \
      --include='*_group_*.csv' --exclude='*' \
      "$HOST:$rdir/" "$ldir/" || echo "WARN: missing $rdir"
  done
}

if [[ $# -eq 0 ]]; then
  set -- h1n1
fi
for d in "$@"; do
  pull_split "$d"
done
echo "Done. Remember: do not commit GISAID HA FASTAs (h3n2)."
