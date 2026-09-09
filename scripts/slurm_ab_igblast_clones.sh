#!/bin/bash
#SBATCH --job-name=ab_igblast
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
# omit --mem: PARCC filter requires exact 5632MB/CPU on genoa-std-mem
#SBATCH --time=48:00:00
#SBATCH --array=0-63
#SBATCH --output=logs/ab_igblast_%A_%a.log
#SBATCH --error=logs/ab_igblast_%A_%a.log
#
# SCAFFOLD — antibody IgBLAST + (optional) Change-O clone assign on shards from
# scripts/prepare_ab_fasta.py. Does NOT claim trees are done.
#
# Prerequisites (edit paths below for Betty install):
#   1. data/ab_prep/shards/ab_shard_XXXX.fasta from slurm_ab_prepare_fasta.sh
#   2. IgBLAST binary + IMGT germline DBs (human IG)
#   3. Optional: Change-O (AssignGenes.py / MakeDb.py / DefineClones.py)
#
# Usage:
#   # After prep finishes:
#   sbatch scripts/slurm_ab_igblast_clones.sh
#   # Fewer shards than array size is OK (missing shard → skip).
#   # Override array to match N_SHARDS:
#   sbatch --array=0-3 scripts/slurm_ab_igblast_clones.sh   # smoke
#
# After all array tasks:
#   # Merge TSVs, filter clones size ∈ [16,64], then (NOT in this script):
#   #   scripts/slurm_ab_build_trees.sh  (MAFFT + FastTree) — still to write
#   # Trees → data/ab_clones/{train,val,test}/ — NOT built by this job.

set -euo pipefail

PREP_DIR="${PREP_DIR:-data/ab_prep}"
OUT_DIR="${OUT_DIR:-data/ab_igblast}"
# Edit these on Betty after installing IgBLAST + IMGT refs:
IGBLAST_BIN="${IGBLAST_BIN:-igblastn}"
GERMLINE_DIR="${GERMLINE_DIR:-/path/to/imgt/human/Ig}"   # V/D/J fasta dirs
ORGANISM="${ORGANISM:-human}"
# IgBLAST is typically nucleotide; our prep shards are AA. For AA-only input,
# either (a) back-translate stub / use ANARCI instead, or (b) point INPUT_MODE=aa
# and skip IgBLAST until NT repertoire dumps are wired. Default: dry scaffolding.
INPUT_MODE="${INPUT_MODE:-scaffold}"   # scaffold | skip

SHARD_ID=$(printf "%04d" "${SLURM_ARRAY_TASK_ID:-0}")
SHARD_FASTA="$PREP_DIR/shards/ab_shard_${SHARD_ID}.fasta"
TASK_OUT="$OUT_DIR/shard_${SHARD_ID}"

mkdir -p logs "$TASK_OUT"

cd ~/DiscreteTreeFlows

echo "Start: $(date)  shard=$SHARD_ID  mode=$INPUT_MODE"
echo "Shard FASTA: $SHARD_FASTA"

if [[ ! -f "$SHARD_FASTA" ]]; then
  echo "WARN: missing shard $SHARD_FASTA — nothing to do"
  exit 0
fi

NSEQ=$(grep -c '^>' "$SHARD_FASTA" || true)
echo "n_seqs_in_shard=$NSEQ"

if [[ "$INPUT_MODE" == "scaffold" || "$INPUT_MODE" == "skip" ]]; then
  cat > "$TASK_OUT/README.txt" <<EOF
IgBLAST scaffold only (INPUT_MODE=$INPUT_MODE).
Shard: $SHARD_FASTA ($NSEQ seqs)
IGBLAST_BIN=$IGBLAST_BIN
GERMLINE_DIR=$GERMLINE_DIR

Next real steps (manual / follow-up scripts):
  1. Install IgBLAST + IMGT human IG germlines; set IGBLAST_BIN / GERMLINE_DIR.
  2. If shards are amino-acid (Homo_sapiens.fasta), either:
       - run ANARCI / AbNumber for V/J + CDR annotation on AA, OR
       - obtain paired NT repertoire and re-run prepare_ab_fasta on NT.
  3. Example IgBLAST (NT) once refs exist:
       \$IGBLAST_BIN -germline_db_V \$GERMLINE_DIR/V \\
         -germline_db_D \$GERMLINE_DIR/D -germline_db_J \$GERMLINE_DIR/J \\
         -organism \$ORGANISM -domain_system imgt -ig_seqtype Ig \\
         -outfmt 19 -query \$SHARD_FASTA -out \$TASK_OUT/igblast.tsv
  4. Change-O: MakeDb.py igblast → DefineClones.py → filter n_leaves∈[16,64].
  5. Tree build: scripts/slurm_ab_build_trees.sh (not written yet).
  6. Do NOT claim Table 6 numbers until Newick + eval exist.

Status: NO IgBLAST run performed by this job.
EOF
  echo "Wrote scaffold README → $TASK_OUT/README.txt"
  echo "Done (scaffold): $(date)"
  exit 0
fi

# Real IgBLAST path (enable with INPUT_MODE=nt once binaries + NT FASTA ready)
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH

if ! command -v "$IGBLAST_BIN" >/dev/null 2>&1; then
  echo "ERROR: IGBLAST_BIN=$IGBLAST_BIN not found on PATH"
  exit 1
fi
if [[ ! -d "$GERMLINE_DIR" ]]; then
  echo "ERROR: GERMLINE_DIR=$GERMLINE_DIR missing — set IMGT refs"
  exit 1
fi

"$IGBLAST_BIN" \
  -germline_db_V "$GERMLINE_DIR/V" \
  -germline_db_D "$GERMLINE_DIR/D" \
  -germline_db_J "$GERMLINE_DIR/J" \
  -organism "$ORGANISM" \
  -domain_system imgt \
  -ig_seqtype Ig \
  -outfmt 19 \
  -num_threads "${SLURM_CPUS_PER_TASK:-8}" \
  -query "$SHARD_FASTA" \
  -out "$TASK_OUT/igblast.tsv"

echo "IgBLAST done: $(date) → $TASK_OUT/igblast.tsv"
echo "Next: Change-O MakeDb/DefineClones + slurm_ab_build_trees (not yet)."
