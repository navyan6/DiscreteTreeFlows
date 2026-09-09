#!/usr/bin/env bash
# Local → Betty: sync lean antibody_benchmark code + frozen trees, then submit full run.
# Requires working SSH (kinit + Duo mux). NO git push.
set -euo pipefail
REPO="${REPO:-/Users/navyanori/Documents/GitHub/DiscreteTreeFlows}"
HOST="${HOST:-nnori@login.betty.parcc.upenn.edu}"
REMOTE_ROOT='~/DiscreteTreeFlows'
SSH=(ssh -o ControlMaster=auto -o ControlPath=~/.ssh/cm-%r@%h:%p -o ControlPersist=30m
     -o BatchMode=yes -o ConnectTimeout=30 "$HOST")

cd "$REPO"
chmod +x antibody_benchmark/scripts/slurm_ab_full_rollout.sh \
  antibody_benchmark/scripts/betty_submit_full.sh \
  .betty_ab_full_submit.sh

echo "=== rsync lean code + frozen trees ==="
rsync -avh --relative \
  --exclude '__pycache__/' \
  --exclude '*.pyc' \
  --exclude '.DS_Store' \
  --exclude 'antibody_benchmark/data/raw/' \
  --exclude 'antibody_benchmark/results/samples/' \
  ./antibody_benchmark/configs \
  ./antibody_benchmark/scripts \
  ./antibody_benchmark/models \
  ./antibody_benchmark/metrics \
  ./antibody_benchmark/rollout \
  ./antibody_benchmark/tests \
  ./antibody_benchmark/codon.py \
  ./antibody_benchmark/trees.py \
  ./antibody_benchmark/__init__.py \
  ./antibody_benchmark/README.md \
  ./antibody_benchmark/DATA_AUDIT.md \
  ./antibody_benchmark/data/processed \
  ./antibody_benchmark/results/smoke_reports \
  "$HOST:$REMOTE_ROOT/"

# Ensure scripts executable remotely + preflight + submit
"${SSH[@]}" "bash -lc '
set -euo pipefail
if type module >/dev/null 2>&1; then module load slurm 2>/dev/null || true; fi
export PATH=/cm/local/apps/slurm/current/bin:\$PATH
cd \$HOME/DiscreteTreeFlows
chmod +x antibody_benchmark/scripts/slurm_ab_full_rollout.sh antibody_benchmark/scripts/betty_submit_full.sh
# Ensure cosine ckpt symlink
LABHOME=/vast/projects/pranam/lab/nnori
mkdir -p antibody_benchmark/data/raw/cosine/checkpoints
if [[ ! -e antibody_benchmark/data/raw/cosine/checkpoints/cosine_dasm.ckpt ]]; then
  ln -sfn \$LABHOME/antibody_benchmark/cosine_ckpts/cosine_dasm.ckpt \
    antibody_benchmark/data/raw/cosine/checkpoints/cosine_dasm.ckpt
fi
ls -la antibody_benchmark/data/processed/benchmark_trees.jsonl checkpoints/best.pt \
  antibody_benchmark/data/raw/cosine/checkpoints/cosine_dasm.ckpt
bash antibody_benchmark/scripts/betty_submit_full.sh
'"

echo "=== pull submit ids + STATUS ==="
scp -o ControlPath=~/.ssh/cm-%r@%h:%p \
  "$HOST:/vast/projects/pranam/lab/nnori/antibody_benchmark/full_submit_ids.json" \
  /tmp/ab_full_submit_ids.json || true
cat /tmp/ab_full_submit_ids.json 2>/dev/null || true
