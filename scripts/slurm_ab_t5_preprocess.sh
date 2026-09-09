#!/bin/bash
#SBATCH --job-name=ab_t5_prep
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=16
# omit --mem: PARCC filter requires exact 5632MB/CPU on genoa-std-mem (~88G @ 16c)
#SBATCH --time=12:00:00
#SBATCH --output=logs/ab_t5_prep_%j.log
#SBATCH --error=logs/ab_t5_prep_%j.log
#
# Table 5 Ab preprocess (<12h budget): subsample ≤500K → ANARCI → clones →
# FastTree/midpoint root → clonal-lineage holdout split.
#
# Prereqs on Betty:
#   data/Homo_sapiens.fasta  OR  data/ab_t5_500k/sequences.fasta (pre-synced subsample)
#   conda env treesbm with python; ANARCI+HMMER; mafft; FastTree; dendropy
#
# Usage:
#   cd ~/DiscreteTreeFlows && mkdir -p logs
#   # If only full FASTA present:
#   sbatch scripts/slurm_ab_t5_preprocess.sh
#   # Skip subsample if already synced:
#   SKIP_SUBSAMPLE=1 sbatch scripts/slurm_ab_t5_preprocess.sh
#   # Smoke (5K is too sparse for min_size>=8 — use MIN_SIZE=2):
#   MAX_SEQS=5000 MAX_CLONES=20 MIN_SIZE=2 OUT=data/ab_t5_smoke \
#     sbatch scripts/slurm_ab_t5_preprocess.sh

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

MAX_SEQS="${MAX_SEQS:-500000}"
MAX_CLONES="${MAX_CLONES:-400}"
OUT="${OUT:-data/ab_t5_500k}"
SKIP_SUBSAMPLE="${SKIP_SUBSAMPLE:-0}"
SKIP_ANARCI="${SKIP_ANARCI:-0}"
INPUT="${INPUT:-data/Homo_sapiens.fasta}"
MIN_SIZE="${MIN_SIZE:-16}"
MAX_SIZE="${MAX_SIZE:-64}"
ANARCI_BATCH="${ANARCI_BATCH:-500}"
ANARCI_JOBS="${ANARCI_JOBS:-8}"

echo "Start: $(date) host=$(hostname)"
echo "MAX_SEQS=$MAX_SEQS OUT=$OUT SKIP_SUBSAMPLE=$SKIP_SUBSAMPLE SKIP_ANARCI=$SKIP_ANARCI"
echo "MIN_SIZE=$MIN_SIZE MAX_SIZE=$MAX_SIZE ANARCI_BATCH=$ANARCI_BATCH ANARCI_JOBS=$ANARCI_JOBS"

if [[ "$SKIP_SUBSAMPLE" != "1" ]]; then
  if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: missing $INPUT — rsync Homo_sapiens.fasta or prebuilt $OUT/sequences.fasta"
    exit 1
  fi
  $PYTHON scripts/prepare_ab_t5_subsample.py \
    --input "$INPUT" \
    --out-dir "$OUT" \
    --max-seqs "$MAX_SEQS" \
    --chains heavy \
    --prefer-paired
else
  echo "Skipping subsample; expecting $OUT/sequences.fasta"
  [[ -f "$OUT/sequences.fasta" ]] || { echo "ERROR: missing $OUT/sequences.fasta"; exit 1; }
fi

if [[ "$SKIP_ANARCI" != "1" ]]; then
  if ! $PYTHON -c "import anarci" 2>/dev/null \
     && ! command -v ANARCI >/dev/null 2>&1 \
     && ! command -v anarci >/dev/null 2>&1; then
    echo "ERROR: ANARCI not in treesbm env. Install then re-run with SKIP_SUBSAMPLE=1"
    echo "  conda install -y -n treesbm -c bioconda -c conda-forge anarci hmmer"
    exit 2
  fi
  $PYTHON scripts/ab_anarci_annotate.py \
    --fasta "$OUT/sequences.fasta" \
    --out "$OUT/anarci.tsv" \
    --batch-size "$ANARCI_BATCH" \
    --n-jobs "$ANARCI_JOBS"
else
  echo "SKIP_ANARCI=1 — expecting $OUT/anarci.tsv"
  [[ -f "$OUT/anarci.tsv" ]] || { echo "ERROR: missing anarci.tsv"; exit 1; }
fi

$PYTHON scripts/ab_define_clones.py \
  --metadata "$OUT/metadata.tsv" \
  --anarci "$OUT/anarci.tsv" \
  --fasta "$OUT/sequences.fasta" \
  --out-dir "$OUT/clones" \
  --min-size "$MIN_SIZE" \
  --max-size "$MAX_SIZE"

N_CLONES=$($PYTHON - <<PY
import json
from pathlib import Path
p=Path("$OUT/clones/summary.json")
print(json.loads(p.read_text()).get("n_kept_clones", 0) if p.is_file() else 0)
PY
)
echo "n_kept_clones=$N_CLONES"
if [[ "$N_CLONES" -eq 0 ]]; then
  echo "WARN: zero clones — stopping before FastTree. See $OUT/clones/summary.json"
  echo "Done (partial): $(date)"
  exit 3
fi

$PYTHON scripts/ab_build_clone_trees.py \
  --clones-dir "$OUT/clones" \
  --out-dir "$OUT/trees_all" \
  --max-clones "$MAX_CLONES"

$PYTHON scripts/ab_family_holdout_split.py \
  --trees-dir "$OUT/trees_all" \
  --clones-index "$OUT/clones/clones_index.tsv" \
  --out-base data/ab_clones \
  --donor-disjoint

echo "Done: $(date)"
echo "Next: baselines via scripts/eval_ab_maturation.py (plm_prior); Neutral SHM / AR still stubs"
echo "  sbatch scripts/slurm_ab_t5_eval.sh"
echo "  $PYTHON scripts/eval_ab_maturation.py --data data/ab_clones/test --methods plm_prior --out results/ab_t5/eval.json"
