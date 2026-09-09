#!/bin/bash
# One-command Betty submit for covid_v4_mutrec enrichment @ mrs=0.3,0.5.
# Use after training finishes (or anytime best.pt exists).
#
#   bash scripts/run_eval_covid_v4_ready.sh
set -euo pipefail
cd ~/DiscreteTreeFlows
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH

if [[ ! -f checkpoints/covid_v4_mutrec/best.pt ]]; then
  echo "ERROR: checkpoints/covid_v4_mutrec/best.pt missing"
  exit 1
fi
mkdir -p logs
JOB=$(sbatch --qos=mig-max --parsable scripts/slurm_eval_covid_v4.sh)
echo "Submitted eval_covid_v4 -> $JOB"
squeue -j "$JOB" -o "%.18i %.12P %.20j %.8T %.10M %.10L %R %q"
