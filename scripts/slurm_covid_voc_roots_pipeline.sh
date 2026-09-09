#!/bin/bash
#SBATCH --job-name=voc_roots_pipeline
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=24:00:00
#SBATCH --output=logs/voc_roots_pipeline_%j.log
#SBATCH --error=logs/voc_roots_pipeline_%j.log
#
# Build trees for VOC-origin populations whose ROOTS the model has never seen
# (results/voc_threat_panel/LEAKAGE.md). Three stages, all CPU:
#   1. nextclade  -> Spike CDS per population
#   2. grouping   -> temporally stratified trees spanning each emergence window
#   3. run_all_groups (mafft -> fasttree -> augur refine -> ancestral -> translate)
#
# Sequences were pulled locally by download_sars2_voc_roots_ncbi.py and are
# verified disjoint from every accession in data/covid/{train,val,test}.
#
# Prereq (data/ is gitignored, so rsync the pull up first):
#   rsync -avP data/covid_voc_roots/raw/ <cluster>:~/DiscreteTreeFlows/data/covid_voc_roots/raw/
#
# Resumable: run_all_groups skips completed stages, so just re-submit on timeout.

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs
shopt -s nullglob

echo "Start: $(date)"

RAW=data/covid_voc_roots/raw
files=("$RAW"/*_covid_seqs.fasta)
if [ ${#files[@]} -eq 0 ]; then
    echo "No $RAW/*_covid_seqs.fasta -- rsync the pull up first (see header)."
    exit 1
fi

echo "=== stage 1: Spike extraction ==="
for f in "${files[@]}"; do
    echo "--- $f ---"
    $PYTHON scripts/covid_extract_spike.py "$f" --jobs "$SLURM_CPUS_PER_TASK" || {
        echo "FAILED spike extraction: $f"; exit 1; }
done

echo "=== stage 2: grouping into stratified trees ==="
$PYTHON scripts/prepare_covid_voc_roots.py --n-trees 3 --group-size 300 \
    --target-frac 0.4 --require-spike || exit 1

echo "=== stage 3: alignment + tree + ancestral reconstruction ==="
fail=0
for d in data/covid_voc_roots/*/; do
    pop=$(basename "$d")
    [ "$pop" = "raw" ] && continue
    [ -z "$(echo "$d"vocroots*_group_*.fasta)" ] && continue
    echo "--- $pop ---"
    if $PYTHON scripts/run_all_groups.py --data-dir "data/covid_voc_roots/$pop" \
        --prefix "vocroots$pop" --workers 16 --stop-after translate; then
        echo "OK: $pop"
    else
        echo "FAILED: $pop"
        fail=1
    fi
done

echo "Done: $(date)"
if [ "$fail" -eq 1 ]; then
    echo "One or more populations FAILED -- see 'FAILED:' lines above."
    exit 1
fi
echo "Next: scripts/screen_voc_threat_trees.py over data/covid_voc_roots/<pop>"
