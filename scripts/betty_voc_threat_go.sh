#!/bin/bash
# Double-click or: bash scripts/betty_voc_threat_go.sh
# Requires: kinit + Duo push on first SSH.
set -euo pipefail
cd "$(dirname "$0")/.."
export BETTY="${BETTY:-nnori@login.betty.parcc.upenn.edu}"
export REMOTE="${REMOTE:-~/DiscreteTreeFlows}"
export DRY="${DRY:-0}"
export SUBMIT="${SUBMIT:-1}"

echo "Opening SSH (approve Duo if prompted)…"
ssh -o ControlMaster=auto -o ControlPath="$HOME/.ssh/cm-%r@%h:%p" -o ControlPersist=4h \
  "$BETTY" 'echo Betty OK; hostname'

# Ensure COVID test/val symlinks (idempotent)
ssh "$BETTY" 'bash -s' <<'EOF'
LAB=/vast/projects/pranam/lab/nnori/DiscreteTreeFlows_data/covid
CD=~/DiscreteTreeFlows/data/covid
mkdir -p "$CD"
ln -sfn "$LAB/test" "$CD/test"
ln -sfn "$LAB/val" "$CD/val"
ln -sfn "$LAB/train" "$CD/train"
ln -sfn "$LAB/wt.txt" "$CD/wt.txt"
ln -sfn "$LAB/evescape_spike_rbd.pt" "$CD/evescape_spike_rbd.pt"
EOF

# Ship the new unseen-root pull (~140MB, gitignored so it must go over rsync) and
# kick off the tree build. These roots are what the leaked origin cases are being
# replaced with -- see results/voc_threat_panel/LEAKAGE.md.
if [ -d data/covid_voc_roots/raw ]; then
  echo "Syncing data/covid_voc_roots/raw …"
  rsync -avP -e "ssh -o ControlPath=$HOME/.ssh/cm-%r@%h:%p" \
    data/covid_voc_roots/raw/ "$BETTY:$REMOTE/data/covid_voc_roots/raw/"
  rsync -avP -e "ssh -o ControlPath=$HOME/.ssh/cm-%r@%h:%p" \
    scripts/download_sars2_voc_roots_ncbi.py \
    scripts/prepare_covid_voc_roots.py \
    scripts/slurm_covid_voc_roots_pipeline.sh \
    "$BETTY:$REMOTE/scripts/"
  if [ "$SUBMIT" = "1" ]; then
    ssh "$BETTY" "bash -lc \"cd $REMOTE && sbatch scripts/slurm_covid_voc_roots_pipeline.sh\""
  fi
fi

bash scripts/betty_sync_and_submit_voc_threat.sh
echo "Pull submit ids when done:"
echo "  scp \$BETTY:~/DiscreteTreeFlows/results/voc_threat_panel/submit_ids.json results/voc_threat_panel/"
echo "  scp \$BETTY:~/DiscreteTreeFlows/results/voc_threat_panel/baseline_submit_ids.json results/voc_threat_panel/"
