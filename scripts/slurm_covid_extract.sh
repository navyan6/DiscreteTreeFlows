#!/bin/bash
#SBATCH --job-name=covid_extract
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/covid_extract_%j.log
#SBATCH --error=logs/covid_extract_%j.log
#
# Spike-gene extraction (CPU): nextclade-align every raw whole-genome FASTA in
# data/covid/train/*_covid_seqs.fasta to reference coordinates, then slice the
# Spike CDS. Run once per raw file (new files added later just need this
# re-run, existing *_aligned_genome.fasta are skipped if already present).
#
# One-time setup, if not already done:
#   conda install -n treesbm -c bioconda -c conda-forge nextclade
#   nextclade dataset get --name sars-cov-2 --output-dir data/covid/nextclade_dataset
#
# Prereq (data/ is gitignored -- rsync raw files up from local, they don't
# arrive via git pull):
#   rsync -avP data/covid/train/*_covid_seqs.fasta <cluster>:~/DiscreteTreeFlows/data/covid/train/

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs
shopt -s nullglob

echo "Start: $(date)"

files=(data/covid/train/*_covid_seqs.fasta)
if [ ${#files[@]} -eq 0 ]; then
    echo "No data/covid/train/*_covid_seqs.fasta found -- rsync the raw files up first (see header)."
    exit 1
fi

for f in "${files[@]}"; do
    echo "=== $f ==="
    $PYTHON scripts/covid_extract_spike.py "$f" --jobs "$SLURM_CPUS_PER_TASK"
done

echo "Done: $(date)"
