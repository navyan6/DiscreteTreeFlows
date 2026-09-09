#!/bin/bash
set -euo pipefail
HOST=nnori@login.betty.parcc.upenn.edu
CTRL=/Users/navyanori/.ssh/cm-nnori@login.betty.parcc.upenn.edu:22
LOCAL=/Users/navyanori/Documents/GitHub/DiscreteTreeFlows
OUT=/tmp/betty_migmax_ops_out.txt
{
  echo "=== START $(date) ==="
  # Sync scripts
  rsync -av -e "ssh -o ControlPath=$CTRL -o BatchMode=yes" \
    "$LOCAL/scripts/slurm_eval_covid.sh" \
    "$LOCAL/scripts/slurm_h1n1_leafholdout_eval.sh" \
    "$LOCAL/scripts/slurm_baselines.sh" \
    "$LOCAL/scripts/slurm_validate_transformer.sh" \
    "$LOCAL/scripts/slurm_inference_sweep.sh" \
    "$LOCAL/scripts/slurm_artreeformer.sh" \
    "$HOST:~/DiscreteTreeFlows/scripts/"
  ssh -o ControlPath="$CTRL" -o BatchMode=yes "$HOST" 'bash -s' << 'REMOTE'
set -euo pipefail
cd ~/DiscreteTreeFlows
mkdir -p logs
echo "=== BEFORE ==="
squeue -u nnori -o "%.18i %.12P %.20j %.8u %.2t %.10M %.6D %R %q"
echo "=== CANCEL covid_mutrec 7323331 ==="
scancel 7323331 || true
# Leave 7323287 if running; cancel PD priority jobs for mig-max resubmit
echo "=== CANCEL PD evals for resubmit ==="
for j in 7323288 7323289 7323290 7323291; do
  st=$(squeue -j "$j" -h -o %t 2>/dev/null || true)
  echo "job $j state=$st"
  if [[ "$st" == "PD" || "$st" == "CF" ]]; then
    scancel "$j" && echo "cancelled $j"
  elif [[ -z "$st" ]]; then
    echo "job $j already gone"
  else
    echo "skip cancel $j (state=$st)"
  fi
done
# Also cancel 7323287 only if still PD (not if R)
st=$(squeue -j 7323287 -h -o %t 2>/dev/null || true)
echo "job 7323287 state=$st"
if [[ "$st" == "PD" ]]; then
  scancel 7323287 && echo "cancelled 7323287 for mig-max resubmit"
else
  echo "keeping 7323287 (state=$st)"
fi
sleep 2
echo "=== RESUBMIT with mig-max ==="
# h1n1
J_H1=$(sbatch --qos=mig-max --parsable scripts/slurm_h1n1_leafholdout_eval.sh) && echo "h1n1_lh_eval -> $J_H1"
# baselines with h3n2 ckpt
if [[ -f checkpoints/h3n2_v2/best.pt ]]; then
  CKPT=checkpoints/h3n2_v2/best.pt
elif [[ -f "${LABHOME:-/vast/projects/pranam/lab}/checkpoints_backup/h3n2_v2/best.pt" ]]; then
  CKPT=checkpoints/h3n2_v2/best.pt
else
  CKPT=checkpoints/h3n2_v2/best.pt
fi
J_BL=$(sbatch --qos=mig-max --parsable scripts/slurm_baselines.sh "$CKPT") && echo "baselines -> $J_BL (ckpt=$CKPT)"
# optional: inf_sweep
J_IS=$(sbatch --qos=mig-max --parsable scripts/slurm_inference_sweep.sh) && echo "inf_sweep -> $J_IS"
# optional: artreeformer
J_AR=$(sbatch --qos=mig-max --parsable scripts/slurm_artreeformer.sh) && echo "artreeformer -> $J_AR"
# GraphTF validate if not already queued
if squeue -u nnori -h -n transformer_val | grep -q .; then
  echo "transformer_val already in queue:"
  squeue -u nnori -n transformer_val
else
  # prefer h3n2_v2
  if [[ -f checkpoints/h3n2_v2/best.pt ]] || [[ -f "${LABHOME:-}/checkpoints_backup/h3n2_v2/best.pt" ]]; then
    J_TV=$(sbatch --qos=mig-max --parsable scripts/slurm_validate_transformer.sh checkpoints/h3n2_v2/best.pt data/train results/transformer_val_h3n2 566) && echo "transformer_val -> $J_TV"
  else
    J_TV=$(sbatch --qos=mig-max --parsable scripts/slurm_validate_transformer.sh) && echo "transformer_val (defaults) -> $J_TV"
  fi
fi
# Only resubmit eval_covid if we cancelled it
if ! squeue -j 7323287 -h >/dev/null 2>&1; then
  J_EC=$(sbatch --qos=mig-max --parsable scripts/slurm_eval_covid.sh) && echo "eval_covid -> $J_EC"
else
  echo "eval_covid 7323287 still present — not resubmitting"
fi
echo "=== AFTER ==="
squeue -u nnori -o "%.18i %.12P %.20j %.8u %.2t %.10M %.6D %R %q"
REMOTE
  echo "=== END $(date) ==="
} >"$OUT" 2>&1
cat "$OUT"
