#!/bin/bash
# Pull formed training trees (nwk + anc_aa + bl.json) from Betty.
# Prefer lab storage (home data/* is often a symlink into DiscreteTreeFlows_data).
#
# Usage:
#   bash scripts/pull_formed_trees_from_betty.sh h1n1
#   bash scripts/pull_formed_trees_from_betty.sh h1n1 covid hiv_geo
#   SKIP_ANC_AA=1 bash scripts/pull_formed_trees_from_betty.sh h3n2   # GISAID: topology only
set -euo pipefail
HOST="${BETTY_HOST:-nnori@login.betty.parcc.upenn.edu}"
# Lab path first (real files); home path is fallback if lab missing.
REMOTE_LAB="${BETTY_DATA:-/vast/projects/pranam/lab/nnori/DiscreteTreeFlows_data}"
REMOTE_HOME="${BETTY_ROOT:-/vast/home/n/nnori/DiscreteTreeFlows}/data"
LOCAL="${ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
SKIP_ANC_AA="${SKIP_ANC_AA:-0}"

pull_split() {
  local name="$1"
  for split in train val test; do
    local ldir="$LOCAL/data/$name/$split"
    mkdir -p "$ldir"
    echo "=== $name/$split ==="
    local includes=(
      --include='group_*_rooted.nwk'
      --include='group_*_bl.json'
      --include='SPLIT_PROTOCOL.json'
      --include='*_group_*.csv'
    )
    if [[ "$SKIP_ANC_AA" != "1" ]]; then
      includes+=(--include='group_*_anc_aa.fasta')
    fi
    # Try lab, then home (rsync -L follows remote symlinks if needed).
    local ok=0
    for rbase in "$REMOTE_LAB" "$REMOTE_HOME"; do
      local rdir="$rbase/$name/$split"
      if rsync -avzL "${includes[@]}" --exclude='*' \
          "$HOST:$rdir/" "$ldir/" 2>/tmp/betty_pull_err.$$; then
        ok=1
        break
      fi
    done
    if [[ "$ok" -ne 1 ]]; then
      echo "WARN: missing $name/$split"
      cat /tmp/betty_pull_err.$$ 2>/dev/null || true
    fi
    rm -f /tmp/betty_pull_err.$$
  done
  # SPLIT_PROTOCOL at dataset root
  rsync -avzL --include='SPLIT_PROTOCOL.json' --include='SIZE_REPORT.json' --exclude='*' \
    "$HOST:$REMOTE_LAB/$name/" "$LOCAL/data/$name/" 2>/dev/null || true
}

if [[ $# -eq 0 ]]; then
  set -- h1n1
fi
for d in "$@"; do
  pull_split "$d"
done
echo "Done. Do not commit GISAID EPI_ISL anc_aa FASTAs (h3n2 / data/train)."
