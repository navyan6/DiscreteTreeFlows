#!/bin/bash
# Submit held-out evals for finished epidemic-split checkpoints.
#
# Prereq: best.pt under checkpoints/*epidemic* and precomputed test trees.
# Usage:
#   bash scripts/betty_submit_epidemic_evals.sh
#   EVAL_MRS="0.5" EVAL_MAX_TREES=20 bash scripts/betty_submit_epidemic_evals.sh

set -euo pipefail
export PATH=/cm/local/apps/slurm/current/bin:/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
export SLURM_CONF="${SLURM_CONF:-/cm/shared/apps/slurm/etc/slurm/slurm.conf}"

cd "${ROOT:-$HOME/DiscreteTreeFlows}"
mkdir -p logs checkpoints benchmarks/results/tables

EVAL_MRS="${EVAL_MRS:-0.5}"
EVAL_MAX_TREES="${EVAL_MAX_TREES:-20}"
NSTEPS="${NSTEPS:-100}"
IDS=benchmarks/results/tables/epidemic_eval_job_ids.json
echo "{" > "$IDS.tmp"
first=1

submit_one() {
  local name="$1"
  local ckpt="$2"
  local data="$3"
  local max_len="$4"
  local use_eve="$5"   # 1 = EVEscape enrichment path
  local mask="${6:-}"

  if [[ ! -f "$ckpt" ]]; then
    echo "SKIP $name — missing $ckpt"
    return
  fi
  if [[ ! -d "$data" ]]; then
    echo "SKIP $name — missing data $data"
    return
  fi

  local out="checkpoints/eval_enrichment_${name}_mrs${EVAL_MRS}.json"
  if [[ -f "$out" ]]; then
    echo "SKIP $name — exists $out"
    return
  fi

  local script
  script=$(mktemp logs/epidemic_eval_${name}_XXXX.sh)
  cat > "$script" <<EOF
#!/bin/bash
#SBATCH --job-name=eval_${name}
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/eval_epidemic_${name}_%j.log
#SBATCH --error=logs/eval_epidemic_${name}_%j.log
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:\$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=\${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=\${HF_HOME:-\$LABHOME/hf_cache}
cd ~/DiscreteTreeFlows
mkdir -p logs checkpoints
echo "Start: \$(date) $name"
EXTRA=()
EOF

  if [[ "$use_eve" == "1" ]]; then
    cat >> "$script" <<EOF
EVESCAPE=data/covid/evescape_spike_rbd.pt
test -f "\$EVESCAPE"
EXTRA+=(--evescape "\$EVESCAPE")
MASK="${mask}"
if [[ -n "\$MASK" && -f "\$MASK" ]]; then
  EXTRA+=(--lit-hotspot-mask "\$MASK")
fi
\$PYTHON -u scripts/eval_evescape_enrichment.py \\
  --checkpoint "$ckpt" \\
  --data "$data" \\
  --max-seq-len $max_len \\
  --mutation-rate-scale $EVAL_MRS \\
  --n-steps $NSTEPS \\
  --max-trees $EVAL_MAX_TREES \\
  --out "$out" \\
  "\${EXTRA[@]}"
EOF
  else
    cat >> "$script" <<EOF
MASK="${mask}"
if [[ -n "\$MASK" && -f "\$MASK" ]]; then
  EXTRA+=(--lit-hotspot-mask "\$MASK")
fi
\$PYTHON -u scripts/eval_evescape_enrichment.py \\
  --checkpoint "$ckpt" \\
  --data "$data" \\
  --max-seq-len $max_len \\
  --mutation-rate-scale $EVAL_MRS \\
  --n-steps $NSTEPS \\
  --max-trees $EVAL_MAX_TREES \\
  --out "$out" \\
  "\${EXTRA[@]}"
EOF
  fi

  cat >> "$script" <<EOF
echo "Done: \$(date) wrote $out"
EOF

  local jid
  jid=$(sbatch --parsable "$script")
  echo "  $name -> $jid  ($out)"
  if [[ $first -eq 1 ]]; then first=0; else echo "," >> "$IDS.tmp"; fi
  printf '  "%s": %s' "$name" "$jid" >> "$IDS.tmp"
}

echo "Submitting epidemic evals (mrs=$EVAL_MRS max_trees=$EVAL_MAX_TREES)..."

submit_one covid_v6_epidemic_mutrec \
  checkpoints/covid_v6_epidemic_mutrec/best.pt \
  data/covid_epidemic/test 1280 1 \
  results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt

submit_one h3n2_v4_epidemic \
  checkpoints/h3n2_v4_epidemic/best.pt \
  data/h3n2_epidemic/test 566 0

submit_one h3n2_v4_epidemic_lit \
  checkpoints/h3n2_v4_epidemic_lit/best.pt \
  data/h3n2_epidemic/test 566 0 \
  results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt

submit_one h1n1_v3_epidemic_mutrec \
  checkpoints/h1n1_v3_epidemic_mutrec/best.pt \
  data/h1n1_epidemic/test 566 0

submit_one h1n1_v3_epidemic_lit_mutrec \
  checkpoints/h1n1_v3_epidemic_lit_mutrec/best.pt \
  data/h1n1_epidemic/test 566 0 \
  results/flu_mutfreq_vs_lit/mut_hotspot_mask_h1_nmicrobiol_lit.pt

echo "" >> "$IDS.tmp"
echo "}" >> "$IDS.tmp"
mv "$IDS.tmp" "$IDS"
echo "Wrote $IDS"
cat "$IDS"
