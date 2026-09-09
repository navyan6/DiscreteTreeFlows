#!/bin/bash
set -euo pipefail
HOST=nnori@login.betty.parcc.upenn.edu
CTRL=/Users/navyanori/.ssh/cm-nnori@login.betty.parcc.upenn.edu:22
LOCAL=/Users/navyanori/Documents/GitHub/DiscreteTreeFlows
OUT=/tmp/betty_eval_v5_submit_out.txt
{
  echo "=== START $(date) ==="
  rsync -av -e "ssh -o ControlMaster=auto -o ControlPath=$CTRL -o ControlPersist=10m -o BatchMode=yes -o ConnectTimeout=30" \
    "$LOCAL/scripts/slurm_eval_covid_v4.sh" \
    "$HOST:~/DiscreteTreeFlows/scripts/"
  ssh -o ControlMaster=auto -o ControlPath="$CTRL" -o ControlPersist=10m -o BatchMode=yes -o ConnectTimeout=30 "$HOST" 'bash -s' << 'REMOTE'
set -euo pipefail
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH
cd ~/DiscreteTreeFlows
echo "=== CKPT ==="
ls -lah checkpoints/covid_v5_mutrec/best.pt
echo "=== QUEUE ==="
squeue -u nnori -o "%.18i %.12P %.20j %.2t %.10M %R %q" | head -50
echo "=== PD count ==="
squeue -u nnori -h -t PD | wc -l
echo "=== SUBMIT eval v5 mrs 0.3 0.5 (train NOT cancelled) ==="
JOB=$(CKPT_NAME=covid_v5_mutrec EVAL_MRS="0.3 0.5" SKIP_EXISTING=0 \
  sbatch --qos=mig-max --parsable \
  --job-name=eval_covid_v5 \
  --output=logs/eval_covid_v5_%j.log \
  --error=logs/eval_covid_v5_%j.log \
  scripts/slurm_eval_covid_v4.sh)
echo "JOB_ID=$JOB"
squeue -j "$JOB" -o "%.18i %.12P %.20j %.8T %.10M %.10L %R %q"
echo "=== TRAIN 7350224 ==="
squeue -j 7350224 -o "%.18i %.12P %.20j %.2t %.10M %R %q" 2>/dev/null || echo "train 7350224 not in queue"
REMOTE
  echo "=== DONE $(date) ==="
} 2>&1 | tee "$OUT"
