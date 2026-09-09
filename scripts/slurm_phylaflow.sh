#!/bin/bash
#SBATCH --job-name=phylaflow
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/vast/home/n/nnori/DiscreteTreeFlows/logs/phylaflow_%j.log
#SBATCH --error=/vast/home/n/nnori/DiscreteTreeFlows/logs/phylaflow_%j.log
#
# PhylaFlow (BHV flow matching): NOT root-conditioned forward generation.
# Table-2 rows once pools exist (EXTERNAL.md / BLOCKERS.md):
#   phylaflow         — native: PhylaFlow topo+BL (rescaled to H) + shared JTT seq
#   phylaflow_adapted — same pool topo + shared BL adapter + shared JTT seq
# Pools keep branch lengths by default (phylaflow_sample.py --keep-branch-lengths).
#
# Default QOS: mig-max (avoids lab MaxGRESPerAccount on qos=mig).
# Override: sbatch --qos=<qos> scripts/slurm_phylaflow.sh <N>
#
# Clone + full protocol (paste on Betty after Duo):
#   bash scripts/_paste_betty_phylaflow_native.sh
#   bash scripts/_paste_betty_phylaflow_native.sh submit
#
# This job does NOT invent pools. It either:
#   (A) symlinks an existing phylaflow_N{N}.nwk into benchmarks/external_pools/sampled/, or
#   (B) post-processes PhylaFlow --dump-trees JSON into that pool via
#       benchmarks/external_adapters/phylaflow_sample.py (keeps BLs), or
#   (C) prints the cluster protocol and exits 1 if nothing is ready.
#
# One N per job (same walltime rationale as slurm_artreeformer.sh).
#
# Usage:
#   sbatch --qos=mig-max scripts/slurm_phylaflow.sh 16
#   sbatch --qos=mig-max scripts/slurm_phylaflow.sh 32
#   sbatch --qos=mig-max scripts/slurm_phylaflow.sh 64
#
# Prereq (once):
#   export LABHOME=/vast/projects/pranam/lab/nnori
#   git clone https://github.com/yashaektefaie/PhylaFlow $LABHOME/baselines/PhylaFlow
#   # install PhylaFlow deps + set PHYLAFLOW_* roots (see EXTERNAL.md)
#   cp benchmarks/external_adapters/phylaflow_sample.py \
#       $LABHOME/baselines/PhylaFlow/phylaflow_sample.py

set -e
N="${1:?usage: sbatch slurm_phylaflow.sh <N>  (16, 32, or 64)}"
LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
REPO_DIR="${LABHOME}/baselines/PhylaFlow"
DTF="${DTF:-/vast/home/n/nnori/DiscreteTreeFlows}"
POOL_SRC="${DTF}/benchmarks/external_pools"
POOL_OUT="${DTF}/benchmarks/external_pools/sampled"
OUT="${POOL_OUT}/phylaflow_N${N}.nwk"
ADAPTER="${DTF}/benchmarks/external_adapters/phylaflow_sample.py"
mkdir -p "$POOL_OUT" "${DTF}/logs"

export PATH="${PATH:-}"
# Prefer treesbm for the post-processor (ete3); PhylaFlow env also fine if active.
if [ -x /vast/home/n/nnori/.conda/envs/treesbm/bin/python ]; then
    PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    echo "ERROR: no python found"
    exit 1
fi

echo "=== PhylaFlow pool builder (N=${N}) → native phylaflow + optional adapted ==="
echo "NOT root-conditioned forward gen — official sampler dumps → pool."
echo "LABHOME=$LABHOME  REPO_DIR=$REPO_DIR"
echo "QOS note: submit with --qos=mig-max (default in this script header)."

if [ ! -d "$REPO_DIR" ]; then
    echo "ERROR: PhylaFlow clone not found at $REPO_DIR"
    echo "  Clone (Betty paste): bash scripts/_paste_betty_phylaflow_native.sh"
    echo "  Or: git clone https://github.com/yashaektefaie/PhylaFlow $REPO_DIR"
    exit 1
fi

echo "=== Expected benchmark pool ==="
echo "  $OUT"
echo "  (one newick topology per line; bare shape, leaf labels anonymized 0..N-1)"

# (A) already-produced pools somewhere on disk → symlink into sampled/
CANDIDATES=(
    "${REPO_DIR}/samples/treesbm_N${N}/phylaflow_N${N}.nwk"
    "${REPO_DIR}/outputs/phylaflow_N${N}.nwk"
    "${POOL_SRC}/phylaflow_N${N}.nwk"
    "${POOL_OUT}/phylaflow_N${N}.nwk"
)
for src in "${CANDIDATES[@]}"; do
    if [ -f "$src" ]; then
        echo "=== Found existing pool: $src ==="
        if [ "$(realpath "$src")" = "$(realpath -m "$OUT")" ]; then
            echo "Already in place: $OUT  ($(wc -l < "$OUT") lines)"
        else
            ln -sfn "$(realpath "$src")" "$OUT"
            echo "Linked -> $OUT  ($(wc -l < "$OUT") lines)"
        fi
        echo "Done: $(date)"
        exit 0
    fi
done

# (B) PhylaFlow tree dumps → convert with thin adapter (no fake trees)
DUMP_CANDIDATES=(
    "${REPO_DIR}/samples/treesbm_N${N}/tree_dumps"
    "${REPO_DIR}/outputs_h3n2/treesbm_N${N}/tree_dumps"
    "${REPO_DIR}/outputs/treesbm_N${N}/tree_dumps"
    "${POOL_SRC}/phylaflow_dumps_N${N}"
)
for dump in "${DUMP_CANDIDATES[@]}"; do
    if [ -d "$dump" ]; then
        echo "=== Converting PhylaFlow dumps: $dump ==="
        if [ ! -f "$ADAPTER" ]; then
            echo "ERROR: missing adapter $ADAPTER"
            exit 1
        fi
        $PYTHON "$ADAPTER" \
            --input "$dump" \
            --ntips "${N}" \
            --n-samples 300 \
            --keep-branch-lengths \
            --out "$OUT"
        echo "Wrote $OUT  ($(wc -l < "$OUT") lines)"
        echo "Done: $(date)"
        exit 0
    fi
done

# (C) Document protocol — do not invent pools
echo "=== Protocol (manual — H3N2 bank train/sample in PhylaFlow env; mig-max) ==="
cat <<EOF
Fair Table 2: train PhylaFlow on OUR H3N2 train data (same empirical pool as
TreeSBM / n≈97 eval). Do NOT use DS1–8 / ./launch_ds_local.sh ds* for pools.

0. Clone + env (once, on lab volume — not home):
     export LABHOME=${LABHOME}
     mkdir -p \$LABHOME/baselines
     git clone https://github.com/yashaektefaie/PhylaFlow \$LABHOME/baselines/PhylaFlow
     cd \$LABHOME/baselines/PhylaFlow
     pip install -r requirements.txt   # or conda env of your choice
     export PHYLAFLOW_DATA_ROOT=\$LABHOME/baselines/phylaflow_h3n2_data
     export PHYLAFLOW_ARTIFACT_ROOT=\$LABHOME/baselines/phylaflow_h3n2_artifacts
     export PHYLAFLOW_OUTPUT_ROOT=\$LABHOME/baselines/PhylaFlow/outputs_h3n2
     # Or: bash ${DTF}/scripts/_paste_betty_phylaflow_native.sh

1. Export anonymized H3N2 train topologies (treesbm env) — same as ARTree/PhyloVAE:
     cd ${DTF}
     python benchmarks/heldout/export_train_topologies.py \\
         --data-dir data/h3n2/train --N ${N} --per-tree 20 \\
         --out-dir benchmarks/external_pools

2. Build custom PhylaFlow bank at \$PHYLAFLOW_DATA_ROOT/h3n2_N${N}/
   (alignments + posteriors / fixed-path anchors + phyla embeddings for H3N2).
   FORBIDDEN for Table 2: short_run_data_DS1-8 / launch_ds_local.sh ds1…ds8.

3. Train PhylaFlow on that H3N2 bank (own env, GPU, qos=mig-max), size≈${N}:
     cd ${REPO_DIR}
     sbatch --qos=mig-max --gres=gpu:1 --partition=b200-mig45 --wrap \\
       'cd ${REPO_DIR} && python -m run.run --config configs/h3n2_N${N}.yaml'

4. Sample terminal trees with PhylaFlow's own sampler + dump newicks:
     cd ${REPO_DIR}
     python scripts/evaluate_per_dataset_sample_kl.py \\
         --config configs/h3n2_N${N}.yaml \\
         --checkpoint \$PHYLAFLOW_OUTPUT_ROOT/<ckpt>.ckpt \\
         --sample-config <sample_metrics.yaml> \\
         --output-dir samples/treesbm_N${N} \\
         --num-samples 50 --dump-trees

5. Convert dumps → benchmark pool (this job, step B):
     python ${ADAPTER} \\
         --input ${REPO_DIR}/samples/treesbm_N${N}/tree_dumps \\
         --ntips ${N} --n-samples 300 \\
         --out ${OUT}
   Or re-run:
     sbatch --qos=mig-max scripts/slurm_phylaflow.sh ${N}

6. Re-run baselines table (native phylaflow + optional phylaflow_adapted):
     sbatch --qos=mig-max scripts/slurm_baselines.sh checkpoints/h3n2_v2/best.pt

Honesty: still not root-conditioned forward gen; pool = H3N2-trained sampler trees.
EOF

echo ""
echo "No pool or dump found yet — complete H3N2 steps above (do not fake .nwk files)."
exit 1
