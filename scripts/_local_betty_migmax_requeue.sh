#!/bin/bash
set -euo pipefail
HOST=nnori@login.betty.parcc.upenn.edu
LOCAL=/Users/navyanori/Documents/GitHub/DiscreteTreeFlows
SSH=(ssh -o BatchMode=yes -o ConnectTimeout=30 -o ControlMaster=auto -o ControlPath=/Users/navyanori/.ssh/cm-%r@%h:%p -o ControlPersist=10m)
RSYNC_E="ssh -o BatchMode=yes -o ControlMaster=auto -o ControlPath=/Users/navyanori/.ssh/cm-%r@%h:%p"
rsync -av -e "$RSYNC_E" \
  "$LOCAL/scripts/slurm_eval_covid.sh" \
  "$LOCAL/scripts/slurm_h1n1_leafholdout_eval.sh" \
  "$LOCAL/scripts/slurm_baselines.sh" \
  "$LOCAL/scripts/slurm_validate_transformer.sh" \
  "$LOCAL/scripts/slurm_inference_sweep.sh" \
  "$LOCAL/scripts/slurm_artreeformer.sh" \
  "$LOCAL/scripts/_paste_betty_migmax_requeue.sh" \
  "$HOST:~/DiscreteTreeFlows/scripts/"
"${SSH[@]}" "$HOST" 'bash ~/DiscreteTreeFlows/scripts/_paste_betty_migmax_requeue.sh'
