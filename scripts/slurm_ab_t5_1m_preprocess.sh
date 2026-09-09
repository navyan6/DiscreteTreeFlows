#!/bin/bash
#SBATCH --job-name=ab_t5_1m_prep
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=16
# omit --mem: PARCC filter requires exact 5632MB/CPU on genoa-std-mem (~88G @ 16c)
#SBATCH --time=18:00:00
#SBATCH --output=logs/ab_t5_1m_prep_%j.log
#SBATCH --error=logs/ab_t5_1m_prep_%j.log
#
# Antibody redo preprocess: 1M AA sequences from Homo_sapiens.fasta
# (same pipeline as 500k Table-5 path; does NOT touch Rodriguez 82 freeze).
#
# Pipeline: subsample → ANARCI → clones (donor|Vfam|CDR3) → FastTree →
# donor-disjoint clonal-lineage holdout → CDR hotspot mask.
#
# Usage:
#   cd ~/DiscreteTreeFlows && mkdir -p logs
#   sbatch scripts/slurm_ab_t5_1m_preprocess.sh
#
# Env overrides:
#   MAX_SEQS=1000000 OUT=data/ab_t5_1m OUT_BASE=data/ab_clones_1m
#   MAX_CLONES=800 MIN_SIZE=16 MAX_SIZE=64
#   SKIP_SUBSAMPLE=1 SKIP_ANARCI=1

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs results/ab_cdr_mask

MAX_SEQS="${MAX_SEQS:-1000000}"
MAX_CLONES="${MAX_CLONES:-800}"
OUT="${OUT:-data/ab_t5_1m}"
OUT_BASE="${OUT_BASE:-data/ab_clones_1m}"
SKIP_SUBSAMPLE="${SKIP_SUBSAMPLE:-0}"
SKIP_ANARCI="${SKIP_ANARCI:-0}"
INPUT="${INPUT:-data/Homo_sapiens.fasta}"
MIN_SIZE="${MIN_SIZE:-16}"
MAX_SIZE="${MAX_SIZE:-64}"
ANARCI_BATCH="${ANARCI_BATCH:-500}"
ANARCI_JOBS="${ANARCI_JOBS:-8}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-160}"
CDR_OCC="${CDR_OCC:-0.25}"

echo "Start: $(date -Is) host=$(hostname) job=${SLURM_JOB_ID:-local}"
echo "MAX_SEQS=$MAX_SEQS OUT=$OUT OUT_BASE=$OUT_BASE MAX_CLONES=$MAX_CLONES"
echo "MIN_SIZE=$MIN_SIZE MAX_SIZE=$MAX_SIZE ANARCI_BATCH=$ANARCI_BATCH ANARCI_JOBS=$ANARCI_JOBS"
echo "holdout=clonal_lineage donor-disjoint (NOT Rodriguez 82)"

if [[ "$SKIP_SUBSAMPLE" != "1" ]]; then
  if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: missing $INPUT"
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
    echo "ERROR: ANARCI not in treesbm env"
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
  echo "WARN: zero clones — stopping before FastTree"
  exit 3
fi

$PYTHON scripts/ab_build_clone_trees.py \
  --clones-dir "$OUT/clones" \
  --out-dir "$OUT/trees_all" \
  --max-clones "$MAX_CLONES"

$PYTHON scripts/ab_family_holdout_split.py \
  --trees-dir "$OUT/trees_all" \
  --clones-index "$OUT/clones/clones_index.tsv" \
  --out-base "$OUT_BASE" \
  --donor-disjoint

# CDR-weighted hotspot mask (empirical from ANARCI CDRs; stub fallback)
MASK_OUT="results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt"
if [[ -f "$OUT/anarci.tsv" && -f "$OUT/sequences.fasta" ]]; then
  $PYTHON scripts/build_ab_cdr_hotspot_mask.py \
    --fasta "$OUT/sequences.fasta" \
    --anarci "$OUT/anarci.tsv" \
    --max-seq-len "$MAX_SEQ_LEN" \
    --min-occupancy "$CDR_OCC" \
    --out "$MASK_OUT"
else
  $PYTHON scripts/build_ab_cdr_hotspot_mask.py \
    --stub --max-seq-len "$MAX_SEQ_LEN" \
    --out "$MASK_OUT"
fi

# Pipeline manifest for train/eval handoff
$PYTHON - <<PY
import json
from pathlib import Path
from datetime import datetime, timezone
manifest = {
    "created": datetime.now(timezone.utc).isoformat(),
    "job_id": "${SLURM_JOB_ID:-local}",
    "input": "$INPUT",
    "max_seqs": int("$MAX_SEQS"),
    "out": "$OUT",
    "out_base": "$OUT_BASE",
    "holdout": "clonal_lineage donor-disjoint (ab_family_holdout_split)",
    "holdout_note": "NOT Rodriguez 82; antibody_benchmark freeze stays separate for CoSiNE/DASM table",
    "cdr_mask": "$MASK_OUT",
    "max_seq_len": int("$MAX_SEQ_LEN"),
    "n_kept_clones": int("$N_CLONES"),
    "max_clones_built": int("$MAX_CLONES"),
    "asr_note": "majority_leaf_vote_placeholder — germline ASR not wired",
    "next": "sbatch --dependency=afterok:${SLURM_JOB_ID:-0} scripts/slurm_ab_oas_treesbm_train.sh",
}
Path("results/ab_cdr_mask/ab_t5_1m_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps(manifest, indent=2))
PY

echo "Done: $(date -Is)"
echo "Trees: $OUT_BASE/{train,val,test}"
echo "CDR mask: $MASK_OUT"
echo "Next: sbatch --dependency=afterok:${SLURM_JOB_ID} scripts/slurm_ab_oas_treesbm_train.sh"
