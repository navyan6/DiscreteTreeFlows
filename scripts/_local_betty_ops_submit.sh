#!/bin/bash
# From laptop (needs live Kerberos + Duo ControlMaster, or open SSH already):
#   bash scripts/_local_betty_ops_submit.sh
# Syncs track_a / siteaa / GraphTF J.1 scripts and runs paste submitters on Betty.
set -euo pipefail
HOST=nnori@login.betty.parcc.upenn.edu
LOCAL=/Users/navyanori/Documents/GitHub/DiscreteTreeFlows
CTRL=/Users/navyanori/.ssh/cm-%r@%h:%p
SSH=(ssh -o BatchMode=yes -o ConnectTimeout=30 -o ControlMaster=auto -o ControlPath="$CTRL" -o ControlPersist=30m)
RSYNC_E="ssh -o BatchMode=yes -o ControlMaster=auto -o ControlPath=$CTRL -o ControlPersist=30m"

rsync -av -e "$RSYNC_E" \
  "$LOCAL/scripts/slurm_track_a.sh" \
  "$LOCAL/scripts/_paste_betty_track_a.sh" \
  "$LOCAL/scripts/_paste_betty_siteaa_reeval.sh" \
  "$LOCAL/scripts/_paste_betty_graphtf_j1.sh" \
  "$LOCAL/scripts/slurm_eval_covid_v4.sh" \
  "$LOCAL/scripts/slurm_validate_transformer.sh" \
  "$LOCAL/scripts/eval_evescape_enrichment.py" \
  "$LOCAL/scripts/eval_single_tree.py" \
  "$HOST:~/DiscreteTreeFlows/scripts/"

rsync -av -e "$RSYNC_E" \
  "$LOCAL/benchmarks/track_a.py" \
  "$LOCAL/benchmarks/generation.py" \
  "$HOST:~/DiscreteTreeFlows/benchmarks/"

rsync -av -e "$RSYNC_E" \
  "$LOCAL/benchmarks/metrics/sequences.py" \
  "$HOST:~/DiscreteTreeFlows/benchmarks/metrics/"

"${SSH[@]}" "$HOST" 'bash -s' << 'REMOTE'
set -euo pipefail
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH
cd ~/DiscreteTreeFlows
chmod +x scripts/slurm_track_a.sh scripts/_paste_betty_track_a.sh \
  scripts/_paste_betty_siteaa_reeval.sh scripts/_paste_betty_graphtf_j1.sh
echo "==== SITEAA ===="
bash scripts/_paste_betty_siteaa_reeval.sh
echo "==== TRACK_A ===="
bash scripts/_paste_betty_track_a.sh
echo "==== GRAPHTF J1 ===="
bash scripts/_paste_betty_graphtf_j1.sh
REMOTE
