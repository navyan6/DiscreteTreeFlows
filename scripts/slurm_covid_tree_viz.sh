#!/bin/bash
#SBATCH --job-name=covid_tree_viz
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/covid_tree_viz_%j.log
#SBATCH --error=logs/covid_tree_viz_%j.log
#
# Generate 1–3 Newick trees from a held-out COVID test root for viz
# (generated vs observed).
#
# Fast (default): max_leaves=48 / n_steps=40
# Matched leaf count example:
#   GROUP=40 N_SAMPLES=1 MRS=0.5 MAX_LEAVES=207 N_STEPS=160 GEN_TAG=matched \
#     sbatch scripts/slurm_covid_tree_viz.sh
#
# Usage:
#   sbatch scripts/slurm_covid_tree_viz.sh
#   GROUP=40 N_SAMPLES=3 MRS=0.5 MAX_LEAVES=48 N_STEPS=40 \
#     sbatch scripts/slurm_covid_tree_viz.sh

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
export TRANSFORMERS_CACHE=${TRANSFORMERS_CACHE:-$HF_HOME}
export TORCH_HOME=${TORCH_HOME:-$LABHOME/torch_cache}

cd ~/DiscreteTreeFlows
mkdir -p logs results/covid_tree_viz

CKPT_NAME="${CKPT_NAME:-covid_v7_pmc_hotspot}"
CKPT="checkpoints/${CKPT_NAME}/best.pt"
FALLBACK="checkpoints/covid_v5_mutrec/best.pt"
DATA="${DATA:-data/covid/test}"
GROUP="${GROUP:-40}"          # smallest held-out test group (~207 leaves)
MRS="${MRS:-0.5}"
MAX_LEAVES="${MAX_LEAVES:-48}"
N_STEPS="${N_STEPS:-40}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-1280}"
N_SAMPLES="${N_SAMPLES:-3}"
BRANCH_SCALE="${BRANCH_SCALE:-6.0}"
OUT_DIR="${OUT_DIR:-results/covid_tree_viz}"
# Optional tag → group_XXX_generated_<tag>.nwk (avoids overwriting the fast 48-leaf set)
GEN_TAG="${GEN_TAG:-}"
# Random tip subsample of observed tree (same root); 0 = skip
SUBSAMPLE_N="${SUBSAMPLE_N:-48}"
SUBSAMPLE_SEED="${SUBSAMPLE_SEED:-42}"
REWRITE_README="${REWRITE_README:-1}"

if [[ ! -f "$CKPT" ]]; then
  echo "WARN: $CKPT missing — falling back to $FALLBACK"
  CKPT="$FALLBACK"
  CKPT_NAME="covid_v5_mutrec"
fi
if [[ ! -f "$CKPT" ]]; then
  echo "ERROR: no usable checkpoint"
  exit 1
fi

OBS_NWK="$DATA/group_$(printf '%03d' "$GROUP")_rooted.nwk"
OBS_FA="$DATA/group_$(printf '%03d' "$GROUP")_anc_aa.fasta"
if [[ ! -f "$OBS_NWK" ]]; then
  echo "ERROR: missing observed tree $OBS_NWK"
  exit 1
fi

G3=$(printf '%03d' "$GROUP")
cp -f "$OBS_NWK" "$OUT_DIR/group_${G3}_observed.nwk"
if [[ -f "$OBS_FA" ]]; then
  cp -f "$OBS_FA" "$OUT_DIR/group_${G3}_observed_anc_aa.fasta"
fi

# Size-matched subsample of observed tips (for pairing with fast max_leaves gens)
if [[ "${SUBSAMPLE_N}" -gt 0 ]]; then
  $PYTHON - <<PY
from pathlib import Path
from ete3 import Tree
import random
src = Path("$OUT_DIR/group_${G3}_observed.nwk")
out = Path("$OUT_DIR/group_${G3}_observed_subsample.nwk")
n_keep = int("$SUBSAMPLE_N")
seed = int("$SUBSAMPLE_SEED")
t = Tree(str(src), format=1)
leaves = t.get_leaf_names()
if len(leaves) < n_keep:
    raise SystemExit(f"observed has {len(leaves)} tips < SUBSAMPLE_N={n_keep}")
keep = sorted(random.Random(seed).sample(leaves, n_keep))
t.prune(keep, preserve_branch_length=True)
out.write_text(t.write(format=1) + "\n")
print(f"subsample: {len(leaves)} -> {n_keep} tips (seed={seed}) -> {out}")
PY
fi

echo "Start: $(date)"
echo "ckpt=$CKPT  group=$GROUP  mrs=$MRS  max_leaves=$MAX_LEAVES  n_steps=$N_STEPS  samples=$N_SAMPLES  gen_tag=${GEN_TAG:-none}"

# Read epoch from checkpoint for README
EPOCH=$($PYTHON - <<PY
import torch
c=torch.load("$CKPT", map_location="cpu", weights_only=False)
print(c.get("epoch", "?"))
PY
)

SEEDS=(42 43 44)
GEN_STEM="group_${G3}_generated"
if [[ -n "$GEN_TAG" ]]; then
  GEN_STEM="group_${G3}_generated_${GEN_TAG}"
fi

for i in $(seq 0 $((N_SAMPLES - 1))); do
  SEED=${SEEDS[$i]:-$((42 + i))}
  echo "=== sample $((i+1))/$N_SAMPLES  seed=$SEED ==="
  $PYTHON scripts/eval_single_tree.py \
      --checkpoint "$CKPT" \
      --data "$DATA" \
      --group "$GROUP" \
      --max-seq-len "$MAX_SEQ_LEN" \
      --max-leaves "$MAX_LEAVES" \
      --mutation-rate-scale "$MRS" \
      --n-steps "$N_STEPS" \
      --branch-rate-scale "$BRANCH_SCALE" \
      --seed "$SEED"

  # eval_single_tree writes checkpoints/gen_group{N}.{nwk,fasta}
  SRC_NWK="checkpoints/gen_group${GROUP}.nwk"
  SRC_FA="checkpoints/gen_group${GROUP}.fasta"
  if [[ ! -f "$SRC_NWK" ]]; then
    echo "ERROR: expected $SRC_NWK after generation"
    exit 1
  fi
  if [[ "$i" -eq 0 ]]; then
    cp -f "$SRC_NWK" "$OUT_DIR/${GEN_STEM}.nwk"
    [[ -f "$SRC_FA" ]] && cp -f "$SRC_FA" "$OUT_DIR/${GEN_STEM}.fasta"
  fi
  if [[ -n "$GEN_TAG" ]]; then
    cp -f "$SRC_NWK" "$OUT_DIR/${GEN_STEM}_s${SEED}.nwk"
    [[ -f "$SRC_FA" ]] && cp -f "$SRC_FA" "$OUT_DIR/${GEN_STEM}_s${SEED}.fasta"
  else
    cp -f "$SRC_NWK" "$OUT_DIR/group_${G3}_generated_s${SEED}.nwk"
    [[ -f "$SRC_FA" ]] && cp -f "$SRC_FA" "$OUT_DIR/group_${G3}_generated_s${SEED}.fasta"
  fi
done

# Leaf counts for README
LEAF_COUNTS=$($PYTHON - <<PY
from pathlib import Path
from ete3 import Tree
out = Path("$OUT_DIR")
rows = []
for p in sorted(out.glob("group_${G3}_*.nwk")):
    n = len(Tree(str(p), format=1).get_leaves())
    rows.append(f"| \`{p.name}\` | {n} |")
print("\n".join(rows))
PY
)

if [[ "$REWRITE_README" == "1" ]]; then
cat > "$OUT_DIR/README.md" <<EOF
# COVID tree viz artifacts

Paired **observed** vs **model-generated** Newick trees for a held-out COVID test group.
Same group / same root for fair comparison; pair by leaf count for side-by-side viz.

| Field | Value |
|-------|-------|
| Group | ${GROUP} (\`group_${G3}\`) |
| Split | \`${DATA}\` (held-out test) |
| Checkpoint | \`${CKPT}\` (\`${CKPT_NAME}\`, epoch ${EPOCH}) |
| Mutation rate scale (mrs) | ${MRS} |
| This job n_steps | ${N_STEPS} |
| This job max_leaves | ${MAX_LEAVES} |
| This job GEN_TAG | ${GEN_TAG:-"(untagged / fast set)"} |
| branch_rate_scale | ${BRANCH_SCALE} |
| max_seq_len | ${MAX_SEQ_LEN} |
| Seeds | ${SEEDS[*]:0:${N_SAMPLES}} |
| Job | SLURM \${SLURM_JOB_ID:-local} @ \$(date -u +%Y-%m-%dT%H:%MZ) |

## Recommended viz pairs

1. **Fast / size-matched (immediate):** \`group_${G3}_observed_subsample.nwk\` ↔ \`group_${G3}_generated.nwk\`
   (both ~48 tips; subsample is random tip prune of observed, seed=${SUBSAMPLE_SEED})
2. **Full-size matched (preferred when available):** \`group_${G3}_observed.nwk\` ↔ \`group_${G3}_generated_matched.nwk\`
   (both ~observed leaf count; regenerate with \`GEN_TAG=matched MAX_LEAVES≈207\`)

## Leaf counts

| File | Leaves |
|------|--------|
${LEAF_COUNTS}

## Files

- \`group_${G3}_observed.nwk\` — full ground-truth rooted tree
- \`group_${G3}_observed_subsample.nwk\` — observed pruned to \`${SUBSAMPLE_N}\` random tips (seed ${SUBSAMPLE_SEED})
- \`group_${G3}_generated.nwk\` — fast primary sample (seed 42; capped max_leaves)
- \`group_${G3}_generated_s{SEED}.nwk\` — additional fast samples
- \`group_${G3}_generated_matched.nwk\` — leaf-count-matched generation (when present)
- optional \`.fasta\` sidecars with node AA sequences

## Quick viz

\`\`\`bash
# Preferred full-size pair (when matched gen exists)
python - <<'PY'
from ete3 import Tree
obs = Tree("group_${G3}_observed.nwk", format=1)
gen = Tree("group_${G3}_generated_matched.nwk", format=1)
print(len(obs), len(gen))
PY

# Fast size-matched pair
python - <<'PY'
from ete3 import Tree
obs = Tree("group_${G3}_observed_subsample.nwk", format=1)
gen = Tree("group_${G3}_generated.nwk", format=1)
print(len(obs), len(gen))
PY

# Or upload both .nwk files to https://itol.embl.de/
\`\`\`
EOF

# Expand SLURM_JOB_ID / date placeholders
$PYTHON - <<PY
from pathlib import Path
import os, datetime
p = Path("$OUT_DIR/README.md")
text = p.read_text()
text = text.replace("\${SLURM_JOB_ID:-local}", os.environ.get("SLURM_JOB_ID", "local"))
text = text.replace("\$(date -u +%Y-%m-%dT%H:%MZ)", datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%MZ"))
p.write_text(text)
PY
fi

echo "Artifacts in $OUT_DIR:"
ls -la "$OUT_DIR"
echo "Done: $(date)"
