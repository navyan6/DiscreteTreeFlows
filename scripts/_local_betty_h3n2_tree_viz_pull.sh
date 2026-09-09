#!/bin/bash
# Wait for Betty h3n2_tree_viz job, then rsync artifacts locally.
# Run from laptop (opens interactive SSH / Duo if needed).
set -euo pipefail
HOST=nnori@login.betty.parcc.upenn.edu
CTRL=/Users/navyanori/.ssh/cm-nnori@login.betty.parcc.upenn.edu:22
LOCAL=/Users/navyanori/Documents/GitHub/DiscreteTreeFlows
JOB="${JOB:-7439927}"
OUT=/tmp/betty_h3n2_tree_viz_pull_out.txt

{
  echo "=== START $(date) ==="
  echo "Ensuring ControlMaster (Duo if needed)..."
  ssh -o ControlMaster=auto -o ControlPath="$CTRL" -o ControlPersist=120m \
    -o ServerAliveInterval=30 -o ServerAliveCountMax=10 \
    -o ConnectTimeout=60 "$HOST" 'echo connected; hostname; date'

  export SSH_OPTS="-o ControlPath=$CTRL -o BatchMode=yes -o ConnectTimeout=30"
  remote() { ssh $SSH_OPTS "$HOST" "$@"; }

  echo "=== Waiting for job $JOB ==="
  for i in $(seq 1 120); do
    st=$(remote 'export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH; export SLURM_CONF=/cm/shared/apps/slurm/etc/slurm/slurm.conf; squeue -j '"$JOB"' -h -o %t 2>/dev/null || true')
    if [[ -z "$st" ]]; then
      echo "job $JOB left queue at $(date)"
      break
    fi
    echo "$(date +%H:%M:%S) state=$st"
    # show progress crumbs
    remote "tail -n 3 ~/DiscreteTreeFlows/logs/h3n2_tree_viz_${JOB}.log 2>/dev/null || true" || true
    sleep 30
  done

  remote 'export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH; export SLURM_CONF=/cm/shared/apps/slurm/etc/slurm/slurm.conf; sacct -j '"$JOB"' -X --format=JobID,JobName%20,State,Elapsed,ExitCode; echo ---; ls -la ~/DiscreteTreeFlows/results/h3n2_tree_viz/'

  mkdir -p "$LOCAL/results/h3n2_tree_viz"
  rsync -av -e "ssh -o ControlPath=$CTRL -o BatchMode=yes" \
    "$HOST:~/DiscreteTreeFlows/results/h3n2_tree_viz/" \
    "$LOCAL/results/h3n2_tree_viz/"

  echo "=== LOCAL FILES ==="
  ls -la "$LOCAL/results/h3n2_tree_viz/"
  echo "=== DONE $(date) ==="
} 2>&1 | tee "$OUT"
