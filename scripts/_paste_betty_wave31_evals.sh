#!/bin/bash
# Wave 3.1: submit h1n1_temporal / covid_temporal / covid_cladeholdout evals.
# Does NOT cancel eval_covid_v5 (7423272). Paste on Betty or: bash scripts/_paste_betty_wave31_evals.sh
set -euo pipefail
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"
export PATH=/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:$PATH
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints

PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
# Betty policy: qos=mig-max; b200-mig45 requires exactly 6 CPUs/GPU
SBATCH_COMMON=(--partition=b200-mig45 --qos=mig-max --gres=gpu:1 --cpus-per-task=6 --mem-per-cpu=8G --time=12:00:00)

submit_one() {
  local name="$1"
  local script="$2"
  if squeue -u nnori -h -n "$name" 2>/dev/null | grep -q .; then
    echo "SKIP: $name already queued/running:"
    squeue -u nnori -n "$name"
    return 0
  fi
  local j
  j=$(sbatch --parsable "${SBATCH_COMMON[@]}" --job-name="$name" \
    --output="logs/${name}_%j.log" --error="logs/${name}_%j.log" "$script")
  echo "JOB_ID_${name}=$j"
}

# --- write ephemeral job bodies ---
cat > scripts/_job_eval_h1n1_temporal_v1.sh << 'JOB'
#!/bin/bash
#SBATCH --job-name=eval_h1n1_temp
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
export TRANSFORMERS_CACHE=${TRANSFORMERS_CACHE:-$HF_HOME}
export TORCH_HOME=${TORCH_HOME:-$LABHOME/torch_cache}
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints
CKPT=checkpoints/h1n1_temporal_v1/best.pt
[[ -f "$CKPT" ]] || { echo "ERROR: missing $CKPT"; exit 1; }
echo "Start: $(date)  ckpt=$CKPT"
# Enrichment-style mut/cons + identity; no EVEscape for H1N1
$PYTHON -u scripts/eval_evescape_enrichment.py \
  --checkpoint "$CKPT" \
  --data data/h1n1_temporal/test \
  --max-seq-len 566 \
  --mutation-rate-scale 0.3 \
  --n-steps 100 \
  --max-trees 20 \
  --out checkpoints/eval_enrichment_h1n1_temporal_v1_mrs0.3.json
echo "Done: $(date)"
JOB

cat > scripts/_job_eval_covid_temporal_v1.sh << 'JOB'
#!/bin/bash
#SBATCH --job-name=eval_covid_temp
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
export TRANSFORMERS_CACHE=${TRANSFORMERS_CACHE:-$HF_HOME}
export TORCH_HOME=${TORCH_HOME:-$LABHOME/torch_cache}
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints
CKPT=checkpoints/covid_temporal_v1/best.pt
EVESCAPE=data/covid/evescape_spike_rbd.pt
[[ -f "$CKPT" ]] || { echo "ERROR: missing $CKPT"; exit 1; }
EV_ARGS=()
if [[ -f "$EVESCAPE" ]]; then
  EV_ARGS=(--evescape "$EVESCAPE")
  echo "EVEscape OK: $EVESCAPE"
else
  echo "WARN: no EVEscape tensor; mut/cons only"
fi
echo "Start: $(date)  ckpt=$CKPT"
$PYTHON -u scripts/eval_evescape_enrichment.py \
  --checkpoint "$CKPT" \
  --data data/covid_temporal/test \
  --max-seq-len 1280 \
  "${EV_ARGS[@]}" \
  --mutation-rate-scale 0.3 \
  --n-steps 100 \
  --max-trees 20 \
  --out checkpoints/eval_enrichment_covid_temporal_v1_mrs0.3.json
echo "Done: $(date)"
JOB

cat > scripts/_job_eval_covid_cladeholdout_v1.sh << 'JOB'
#!/bin/bash
#SBATCH --job-name=eval_covid_clade
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
export TRANSFORMERS_CACHE=${TRANSFORMERS_CACHE:-$HF_HOME}
export TORCH_HOME=${TORCH_HOME:-$LABHOME/torch_cache}
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints
CKPT=checkpoints/covid_cladeholdout_v1/best.pt
[[ -f "$CKPT" ]] || { echo "ERROR: missing $CKPT"; exit 1; }
echo "Start: $(date)  ckpt=$CKPT"
$PYTHON -u scripts/eval_leaf_holdout.py \
  --data data/covid_cladeholdout \
  --checkpoint "$CKPT" \
  --max-seq-len 1280 \
  --out checkpoints/eval_leaf_holdout_covid_cladeholdout_v1.json
echo "Done: $(date)"
JOB

chmod +x scripts/_job_eval_h1n1_temporal_v1.sh \
         scripts/_job_eval_covid_temporal_v1.sh \
         scripts/_job_eval_covid_cladeholdout_v1.sh

submit_one eval_h1n1_temp   scripts/_job_eval_h1n1_temporal_v1.sh
submit_one eval_covid_temp  scripts/_job_eval_covid_temporal_v1.sh
submit_one eval_covid_clade scripts/_job_eval_covid_cladeholdout_v1.sh

echo "=== squeue ==="
squeue -u nnori -o "%.18i %.9P %.28j %.8u %.2t %.10M %.6D %R %q %.4C"
