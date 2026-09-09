#!/bin/bash
# Paste into active Betty SSH session (cwd ~/DiscreteTreeFlows) OR:
#   bash scripts/_paste_betty_migmax_requeue.sh
set -euo pipefail
cd ~/DiscreteTreeFlows
mkdir -p logs scripts
echo "=== BEFORE ==="
squeue -u nnori -o "%.18i %.12P %.20j %.8u %.2t %.10M %.6D %R %q"

echo "=== CANCEL 7323331 covid_mutrec ==="
scancel 7323331 || true

# Cancel PD jobs for mig-max resubmit; keep running eval_covid
for j in 7323288 7323289 7323290 7323291; do
  st=$(squeue -j "$j" -h -o %t 2>/dev/null || true)
  echo "job $j state=${st:-gone}"
  if [[ "$st" == "PD" || "$st" == "CF" ]]; then
    scancel "$j" && echo "cancelled $j"
  fi
done
st=$(squeue -j 7323287 -h -o %t 2>/dev/null || true)
echo "job 7323287 state=${st:-gone}"
if [[ "$st" == "PD" ]]; then
  scancel 7323287 && echo "cancelled 7323287"
else
  echo "keeping 7323287 (not PD)"
fi
sleep 2

echo "=== SUBMIT mig-max ==="
set +e
J_H1=$(sbatch --qos=mig-max --parsable scripts/slurm_h1n1_leafholdout_eval.sh 2>&1)
echo "h1n1_lh_eval: $J_H1"
J_BL=$(sbatch --qos=mig-max --parsable scripts/slurm_baselines.sh checkpoints/h3n2_v2/best.pt 2>&1)
echo "baselines: $J_BL"
J_IS=$(sbatch --qos=mig-max --parsable scripts/slurm_inference_sweep.sh 2>&1)
echo "inf_sweep: $J_IS"
J_AR=$(sbatch --qos=mig-max --parsable scripts/slurm_artreeformer.sh 2>&1)
echo "artreeformer: $J_AR"
if ! squeue -u nnori -h -n transformer_val 2>/dev/null | grep -q .; then
  J_TV=$(sbatch --qos=mig-max --parsable scripts/slurm_validate_transformer.sh checkpoints/h3n2_v2/best.pt data/h3n2/train results/transformer_val_h3n2_j1 566 2>&1)
  echo "transformer_val: $J_TV"
else
  echo "transformer_val already queued:"
  squeue -u nnori -n transformer_val
fi
if ! squeue -j 7323287 -h >/dev/null 2>&1; then
  J_EC=$(sbatch --qos=mig-max --parsable scripts/slurm_eval_covid.sh 2>&1)
  echo "eval_covid: $J_EC"
fi
set -e
echo "=== AFTER ==="
squeue -u nnori -o "%.18i %.12P %.20j %.8u %.2t %.10M %.6D %R %q"
