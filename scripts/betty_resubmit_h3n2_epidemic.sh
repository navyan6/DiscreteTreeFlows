#!/bin/bash
# Re-submit H3N2 epidemic precompute + train after pipeline/root fix.
set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
cd "${ROOT:-$HOME/DiscreteTreeFlows}"

scancel 8027978 8027981 2>/dev/null || true

j_pipe=$(sbatch --parsable --export=ALL,DATA_ROOT=data/h3n2_epidemic \
  scripts/slurm_epidemic_pipeline.sh)
echo "h3n2 epidemic pipeline (resume) -> $j_pipe"

j_pre=$(sbatch --parsable --dependency=afterok:"$j_pipe" \
  --export=ALL,DATA_ROOT=data/h3n2_epidemic,MAX_SEQ_LEN=566 \
  scripts/slurm_epidemic_precompute.sh)
echo "h3n2 epidemic precompute -> $j_pre"

j1=$(sbatch --parsable --dependency=afterok:"$j_pre" --qos=mig-max \
  --job-name=ep_h3n2_v4_epidemic \
  --export=ALL,DATA_ROOT=data/h3n2_epidemic,CKPT_DIR=checkpoints/h3n2_v4_epidemic_mutrec,MAX_SEQ_LEN=566,HOTSPOT=0 \
  scripts/slurm_epidemic_train_mut_recovery.sh)
echo "  train h3n2_v4_epidemic -> $j1"

j2=$(sbatch --parsable --dependency=afterok:"$j_pre" --qos=mig-max \
  --job-name=ep_h3n2_v4_epidemic_lit \
  --export=ALL,DATA_ROOT=data/h3n2_epidemic,CKPT_DIR=checkpoints/h3n2_v4_epidemic_lit_mutrec,MAX_SEQ_LEN=566,HOTSPOT=1,MUT_HOTSPOT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt \
  scripts/slurm_epidemic_train_mut_recovery.sh)
echo "  train h3n2_v4_epidemic_lit -> $j2"
