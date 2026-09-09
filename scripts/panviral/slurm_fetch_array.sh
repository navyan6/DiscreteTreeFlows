#!/bin/bash
#SBATCH --job-name=pv_fetch
#SBATCH --partition=genoa-std-mem
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --output=logs/panviral/fetch_%A_%a.log

# Stage 2 fan-out: one array task per virus.
#
# Concurrency is capped in the submit line (%N), not here, and it must stay low.
# NCBI rate-limits per source IP, so every running task shares one 3 req/s
# budget -- running 50 tasks at once does not go 50x faster, it just turns the
# whole array into 429 retries. The per-task --sleep is set so that
# concurrency x (1/sleep) lands under the limit.

set -uo pipefail
REPO="${TREESBM_ROOT:-$HOME/DiscreteTreeFlows}"
cd "$REPO"

WORKLIST=${WORKLIST:-data/panviral/worklist.tsv}
MAX_RECORDS=${MAX_RECORDS:-5000}
MIN_LEAVES=${MIN_LEAVES:-30}
MAX_LEAVES=${MAX_LEAVES:-600}
SLEEP=${SLEEP:-1.0}
PY=${TREESBM_PY:-/vast/home/n/nnori/.conda/envs/treesbm/bin/python}

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$WORKLIST")
if [ -z "$line" ]; then
    echo "no work at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

taxid=$(echo "$line" | cut -f1)
slug=$(echo "$line" | cut -f2)
count=$(echo "$line" | cut -f3)

echo "[$(date -Is)] task ${SLURM_ARRAY_TASK_ID}: $slug (taxid $taxid, $count seqs)"

if [ -f "data/panviral/${slug}/manifest.json" ]; then
    echo "already built, skipping"
    exit 0
fi

$PY -W ignore scripts/panviral/fetch_virus_dataset.py \
    --taxid "$taxid" --name "$slug" \
    --max-records "$MAX_RECORDS" \
    --min-leaves "$MIN_LEAVES" --max-leaves "$MAX_LEAVES" \
    --threads "${SLURM_CPUS_PER_TASK:-4}" --sleep "$SLEEP"
rc=$?

# A virus with no usable data is an expected outcome, not a failure: it should
# not mark the array task failed and trigger a pointless requeue.
if [ $rc -ne 0 ]; then
    echo "[$(date -Is)] $slug produced no dataset (rc=$rc) -- recording and continuing"
    mkdir -p data/panviral/_skipped
    echo "$line" >> "data/panviral/_skipped/rc${rc}.tsv"
fi
echo "[$(date -Is)] done $slug"
exit 0
