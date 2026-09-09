#!/bin/bash
# Paste on Betty after Duo SSH (cwd ~/DiscreteTreeFlows).
# Native PhylaFlow Table-2 setup: clone → H3N2 bank train/sample → pool → baselines.
#
# Locked: same empirical pool as table_empirical_main (data/h3n2/test, max-roots=100,
# N∈{16,32}, n_roots≈97 unique). Row label: phylaflow (native topo+BL + shared JTT seq).
#
# Fairness: train/sample PhylaFlow on OUR H3N2 train data (same distribution as
# TreeSBM / ARTreeFormer / PhyloVAE Table 2). Do NOT use PhylaFlow paper DS1–8
# (./launch_ds_local.sh ds*) for Table-2 pools — those reproduce PhylaFlow's
# own tables, not a fair RF/Q/BW/TE comparison on H3N2 held-out roots.
#
# Usage:
#   bash scripts/_paste_betty_phylaflow_native.sh           # clone + status + print next
#   bash scripts/_paste_betty_phylaflow_native.sh submit    # also sbatch pool jobs + table if ready
set -euo pipefail

LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
DTF="${DTF:-$HOME/DiscreteTreeFlows}"
REPO_DIR="${LABHOME}/baselines/PhylaFlow"
MODE="${1:-status}"

cd "$DTF"
mkdir -p logs "$LABHOME/baselines" \
  "$DTF/benchmarks/external_pools/sampled"

echo "=== Native PhylaFlow Table-2 setup (H3N2 bank — NOT DS1–8) ==="
echo "LABHOME=$LABHOME"
echo "DTF=$DTF"
echo "REPO_DIR=$REPO_DIR"
echo "date=$(date)"

# ── 1) Clone (lab volume — not home; home is 50GB) ──────────────────────────
if [ ! -d "$REPO_DIR/.git" ]; then
  echo "=== Cloning PhylaFlow → $REPO_DIR ==="
  git clone https://github.com/yashaektefaie/PhylaFlow "$REPO_DIR"
  (cd "$REPO_DIR" && git rev-parse HEAD | tee "$DTF/logs/phylaflow_clone_pin.txt")
else
  echo "=== PhylaFlow already cloned ==="
  (cd "$REPO_DIR" && git rev-parse HEAD && git remote -v | head -2)
fi

# ── 2) Adapter into PhylaFlow tree ──────────────────────────────────────────
cp -f "$DTF/benchmarks/external_adapters/phylaflow_sample.py" \
  "$REPO_DIR/phylaflow_sample.py"
echo "Copied phylaflow_sample.py → $REPO_DIR/"

# ── 3) Env / data roots (H3N2 custom bank; DS1–8 is out of scope for Table 2) ─
export PHYLAFLOW_OUTPUT_ROOT="${PHYLAFLOW_OUTPUT_ROOT:-$REPO_DIR/outputs_h3n2}"
export PHYLAFLOW_DATA_ROOT="${PHYLAFLOW_DATA_ROOT:-$LABHOME/baselines/phylaflow_h3n2_data}"
export PHYLAFLOW_ARTIFACT_ROOT="${PHYLAFLOW_ARTIFACT_ROOT:-$LABHOME/baselines/phylaflow_h3n2_artifacts}"
mkdir -p "$PHYLAFLOW_OUTPUT_ROOT" "$PHYLAFLOW_DATA_ROOT" "$PHYLAFLOW_ARTIFACT_ROOT"

echo ""
echo "=== H3N2 bank roots (Table 2 — set / populate these) ==="
echo "  export PHYLAFLOW_DATA_ROOT=$PHYLAFLOW_DATA_ROOT"
echo "  export PHYLAFLOW_ARTIFACT_ROOT=$PHYLAFLOW_ARTIFACT_ROOT"
echo "  export PHYLAFLOW_OUTPUT_ROOT=$PHYLAFLOW_OUTPUT_ROOT"
echo "  Expected layout (custom, not DS1–8):"
echo "    \$PHYLAFLOW_DATA_ROOT/h3n2_N16/   # alignments + posterior/path artifacts"
echo "    \$PHYLAFLOW_DATA_ROOT/h3n2_N32/"
echo "    \$PHYLAFLOW_ARTIFACT_ROOT/phyla_embeddings/   # H3N2 cases"
echo "  DS1–8 (short_run_data_DS1-8/, launch_ds_local.sh ds*) is PhylaFlow-paper"
echo "  only — do NOT use for Table-2 pools."

DATA_OK=0
# Prefer an H3N2 bank; tolerate either N16 or N32 present.
if [ -d "${PHYLAFLOW_DATA_ROOT}/h3n2_N16" ] || [ -d "${PHYLAFLOW_DATA_ROOT}/h3n2_N32" ]; then
  DATA_OK=1
  echo "PHYLAFLOW_DATA_ROOT has H3N2 bank dirs: $PHYLAFLOW_DATA_ROOT"
else
  echo "WARNING: H3N2 bank not built yet under $PHYLAFLOW_DATA_ROOT — cannot train/sample."
fi
# Soft warn if someone pointed at DS1–8 by mistake
if [ -d "${PHYLAFLOW_DATA_ROOT}/short_run_data_DS1-8" ] && [ "$DATA_OK" != "1" ]; then
  echo "NOTE: found short_run_data_DS1-8 under DATA_ROOT — that is PhylaFlow-paper data."
  echo "      For Table 2, build h3n2_N{N}/ instead (see EXTERNAL.md)."
fi

# ── 4) Existing dumps / pools? ──────────────────────────────────────────────
echo ""
echo "=== Pool / dump inventory ==="
for N in 16 32 64; do
  pool="$DTF/benchmarks/external_pools/sampled/phylaflow_N${N}.nwk"
  if [ -f "$pool" ]; then
    echo "  OK pool N=$N  lines=$(wc -l < "$pool")  $pool"
  else
    echo "  MISSING pool N=$N  → $pool"
  fi
done
for N in 16 32 64; do
  for d in \
    "$REPO_DIR/samples/treesbm_N${N}/tree_dumps" \
    "$REPO_DIR/outputs_h3n2/treesbm_N${N}/tree_dumps" \
    "$REPO_DIR/outputs/treesbm_N${N}/tree_dumps" \
    "$DTF/benchmarks/external_pools/phylaflow_dumps_N${N}"
  do
    if [ -d "$d" ]; then
      echo "  dump dir: $d  ($(find "$d" -type f | wc -l) files)"
    fi
  done
done

# ── 5) Optional submit ──────────────────────────────────────────────────────
if [ "$MODE" != "submit" ]; then
  cat <<EOF

=== Next steps (automated Table-2 path — NOT DS1–8) ===

A) Setup env + H3N2 banks N=16/32 (mig-max):
     sbatch --qos=mig-max $DTF/scripts/slurm_phylaflow_setup_h3n2.sh
     # → venv at \$LABHOME/baselines/conda_envs/phylaflow
     # → banks at \$PHYLAFLOW_DATA_ROOT/h3n2_N{16,32}/
     # → configs at \$REPO_DIR/configs/h3n2_N{16,32}.yaml
     # Uses existing train_topologies_N*.nwk/.trprobs (export already done).

B) Train (after setup succeeds; one N per job):
     sbatch --qos=mig-max --dependency=afterok:<SETUP_JID> \\
       $DTF/scripts/slurm_phylaflow_train_h3n2.sh 16
     sbatch --qos=mig-max --dependency=afterok:<SETUP_JID> \\
       $DTF/scripts/slurm_phylaflow_train_h3n2.sh 32
     # FORBIDDEN: ./launch_ds_local.sh ds1…ds8

C) After train: sample dumps with PhylaFlow evaluate_per_dataset_sample_kl.py
   --dump-trees → then convert:
     sbatch --qos=mig-max $DTF/scripts/slurm_phylaflow.sh 16
     sbatch --qos=mig-max $DTF/scripts/slurm_phylaflow.sh 32

D) Score Table 2 (n≈97, N∈{16,32} primary):
     sbatch --qos=mig-max $DTF/scripts/slurm_baselines.sh checkpoints/h3n2_v2/best.pt

Job IDs (if submitted): cat $DTF/logs/phylaflow_h3n2_job_ids.txt
Re-run: bash scripts/_paste_betty_phylaflow_native.sh submit
EOF
  exit 0
fi

echo ""
echo "=== submit mode ==="
# Convert any ready dumps → pools
JOBS=()
for N in 16 32 64; do
  pool="$DTF/benchmarks/external_pools/sampled/phylaflow_N${N}.nwk"
  if [ -f "$pool" ]; then
    echo "Pool N=$N already present — skip convert job"
    continue
  fi
  jid=$(sbatch --parsable --qos=mig-max "$DTF/scripts/slurm_phylaflow.sh" "$N" || true)
  if [ -n "${jid:-}" ]; then
    echo "Submitted phylaflow convert N=$N → job $jid"
    JOBS+=("$jid")
  fi
done

# If any pool exists (or convert submitted), queue baselines after convert
have_any=0
for N in 16 32 64; do
  [ -f "$DTF/benchmarks/external_pools/sampled/phylaflow_N${N}.nwk" ] && have_any=1
done

if [ "$have_any" = "1" ]; then
  jid=$(sbatch --parsable --qos=mig-max \
    "$DTF/scripts/slurm_baselines.sh" checkpoints/h3n2_v2/best.pt)
  echo "Submitted baselines (includes native phylaflow row) → job $jid"
elif [ "${#JOBS[@]}" -gt 0 ]; then
  dep=$(IFS=:; echo "${JOBS[*]}")
  jid=$(sbatch --parsable --qos=mig-max --dependency=afterok:${dep} \
    "$DTF/scripts/slurm_baselines.sh" checkpoints/h3n2_v2/best.pt)
  echo "Submitted baselines after convert jobs ($dep) → job $jid"
else
  echo "No pools and convert jobs will exit 1 without dumps."
  echo "Complete H3N2 bank train/sample (steps B–E) first, then re-run submit."
  if [ "$DATA_OK" != "1" ]; then
    echo "BLOCKER: build PHYLAFLOW_DATA_ROOT/h3n2_N{N}/ (not DS1–8)."
  fi
  exit 1
fi

echo "Done submit. Check: squeue -u \$USER | grep -E 'phylaflow|baselines'"
