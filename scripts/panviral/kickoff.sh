#!/usr/bin/env bash
# One entrypoint for the pan-viral data pipeline on Betty (or any SLURM host).
#
#   bash scripts/panviral/kickoff.sh
#
# Submits stage 1 (inventory) → stage 2 (fetch + CDS extract) → stage 3
# (mafft / FastTree / augur / translate). Later stages wait on earlier ones via
# SLURM dependencies, so this is fire-and-forget.
#
# Useful knobs (all optional):
#   TREESBM_ROOT=/path/to/repo          # default: $HOME/DiscreteTreeFlows
#   TREESBM_PY=/path/to/python          # default: treesbm conda on Betty
#   MIN_COUNT=150                       # genomes required to qualify a virus
#   CHAIN_STAGE3=1                      # set 0 to stop after the CDS pull
#   NCBI_API_KEY_FILE=~/.ncbi_api_key   # presence raises NCBI rate limits
#
# Resume is free: inventory caches counts, stage 2 skips viruses with a
# manifest, stage 3 skips splits that already have rooted trees.
set -euo pipefail

REPO="${TREESBM_ROOT:-$HOME/DiscreteTreeFlows}"
cd "$REPO"
mkdir -p logs/panviral data/panviral

if ! command -v sbatch >/dev/null 2>&1; then
    # Betty login shells sometimes lack SLURM on PATH.
    export PATH="/vast/parcc/sw/slurm/bin:${PATH}"
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch not found; run this on a SLURM login node" >&2
    exit 1
fi

export TREESBM_ROOT="$REPO"
export TREESBM_PY="${TREESBM_PY:-/vast/home/n/nnori/.conda/envs/treesbm/bin/python}"
export MIN_COUNT="${MIN_COUNT:-150}"
export CHAIN_STAGE3="${CHAIN_STAGE3:-1}"

echo "repo=$REPO"
echo "python=$TREESBM_PY"
echo "min_count=$MIN_COUNT  chain_stage3=$CHAIN_STAGE3"

inv=$(sbatch --parsable --export=ALL \
      scripts/panviral/slurm_inventory.sh)
echo "stage 1 inventory: $inv"
echo "$inv" > data/panviral/stage1_jobid.txt

s2=$(sbatch --parsable --dependency=afterok:"$inv" --export=ALL \
      scripts/panviral/slurm_stage2_launcher.sh)
echo "stage 2 launcher:  $s2  (waits on $inv)"
echo "$s2" > data/panviral/stage2_launcher_jobid.txt

echo
echo "queued. track with:  squeue -u \$USER"
echo "logs under:          logs/panviral/"
if [ "$CHAIN_STAGE3" = "1" ]; then
    echo "stage 3 submits itself once the stage-2 fetch array finishes."
fi
