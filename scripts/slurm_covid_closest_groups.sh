#!/bin/bash
#SBATCH --job-name=covid_closest
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/covid_closest_%j.log
#SBATCH --error=logs/covid_closest_%j.log
#
# Screen several held-out COVID test groups, rank gen↔obs closeness, then
# regenerate size-matched trees for the top winners (alongside group 40).
#
# Env overrides:
#   CAND_GROUPS="21 8 9 ..."   SCREEN_MAX=96  SCREEN_STEPS=60  SCREEN_SAMPLES=2
#   MATCH_CAP=250  MATCH_STEPS=160  MATCH_SAMPLES=2  PROMOTE_TOP=3
#   MRS=0.5  CKPT_NAME=covid_v7_pmc_hotspot
#
# Note: do NOT use the name GROUPS — bash reserves it (unix gid list).
#
# Usage:
#   sbatch scripts/slurm_covid_closest_groups.sh

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
export TRANSFORMERS_CACHE=${TRANSFORMERS_CACHE:-$HF_HOME}
export TORCH_HOME=${TORCH_HOME:-$LABHOME/torch_cache}

cd ~/DiscreteTreeFlows
mkdir -p logs results/covid_tree_viz results/covid_tree_viz_screen

CKPT_NAME="${CKPT_NAME:-covid_v7_pmc_hotspot}"
CKPT="checkpoints/${CKPT_NAME}/best.pt"
DATA="${DATA:-data/covid/test}"
MRS="${MRS:-0.5}"
BRANCH_SCALE="${BRANCH_SCALE:-6.0}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-1280}"

# Diverse held-out groups excluding 40 (already viz'd). Sizes are ~253–300;
# screen at SCREEN_MAX, then matched at min(obs, MATCH_CAP).
CAND_GROUPS="${CAND_GROUPS:-21 8 9 36 20 1 11 3 15 39}"
SCREEN_MAX="${SCREEN_MAX:-96}"
SCREEN_STEPS="${SCREEN_STEPS:-60}"
SCREEN_SAMPLES="${SCREEN_SAMPLES:-2}"
MATCH_CAP="${MATCH_CAP:-250}"
MATCH_STEPS="${MATCH_STEPS:-160}"
MATCH_SAMPLES="${MATCH_SAMPLES:-2}"
PROMOTE_TOP="${PROMOTE_TOP:-3}"
SKIP_SCREEN="${SKIP_SCREEN:-0}"
SKIP_MATCH="${SKIP_MATCH:-0}"

SCREEN_DIR="${SCREEN_DIR:-results/covid_tree_viz_screen}"
OUT_DIR="${OUT_DIR:-results/covid_tree_viz}"
SEEDS=(42 43 44)

if [[ ! -f "$CKPT" ]]; then
  echo "ERROR: missing checkpoint $CKPT"
  exit 1
fi

obs_leaves() {
  local g="$1"
  $PYTHON - <<PY
from ete3 import Tree
from pathlib import Path
p = Path("$DATA") / f"group_{int('$g'):03d}_rooted.nwk"
print(len(Tree(str(p), format=1).get_leaves()))
PY
}

generate_group() {
  local GROUP="$1"
  local MAX_LEAVES="$2"
  local N_STEPS="$3"
  local N_SAMPLES="$4"
  local GEN_TAG="$5"
  local DEST_DIR="$6"

  local G3
  G3=$(printf '%03d' "$GROUP")
  mkdir -p "$DEST_DIR"

  # Always keep a copy of the full observed tree in DEST_DIR
  cp -f "$DATA/group_${G3}_rooted.nwk" "$DEST_DIR/group_${G3}_observed.nwk"
  if [[ -f "$DATA/group_${G3}_anc_aa.fasta" ]]; then
    cp -f "$DATA/group_${G3}_anc_aa.fasta" "$DEST_DIR/group_${G3}_observed_anc_aa.fasta"
  fi

  local GEN_STEM="group_${G3}_generated"
  if [[ -n "$GEN_TAG" ]]; then
    GEN_STEM="group_${G3}_generated_${GEN_TAG}"
  fi

  echo "=== gen group=$GROUP max_leaves=$MAX_LEAVES n_steps=$N_STEPS samples=$N_SAMPLES tag=${GEN_TAG:-none} ==="
  local i
  for i in $(seq 0 $((N_SAMPLES - 1))); do
    local SEED=${SEEDS[$i]:-$((42 + i))}
    echo "--- sample $((i+1))/$N_SAMPLES seed=$SEED ---"
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

    local SRC_NWK="checkpoints/gen_group${GROUP}.nwk"
    local SRC_FA="checkpoints/gen_group${GROUP}.fasta"
    [[ -f "$SRC_NWK" ]] || { echo "ERROR: missing $SRC_NWK"; exit 1; }

    if [[ "$i" -eq 0 ]]; then
      cp -f "$SRC_NWK" "$DEST_DIR/${GEN_STEM}.nwk"
      [[ -f "$SRC_FA" ]] && cp -f "$SRC_FA" "$DEST_DIR/${GEN_STEM}.fasta"
    fi
    cp -f "$SRC_NWK" "$DEST_DIR/${GEN_STEM}_s${SEED}.nwk"
    [[ -f "$SRC_FA" ]] && cp -f "$SRC_FA" "$DEST_DIR/${GEN_STEM}_s${SEED}.fasta"
  done
}

echo "Start: $(date)"
echo "ckpt=$CKPT mrs=$MRS cand_groups=[$CAND_GROUPS]"
echo "screen: max=$SCREEN_MAX steps=$SCREEN_STEPS samples=$SCREEN_SAMPLES"
echo "match:  cap=$MATCH_CAP steps=$MATCH_STEPS samples=$MATCH_SAMPLES promote=$PROMOTE_TOP"

# ── Phase 1: screen ─────────────────────────────────────────────────────────
if [[ "$SKIP_SCREEN" != "1" ]]; then
  for GROUP in $CAND_GROUPS; do
    N_OBS=$(obs_leaves "$GROUP")
    MAX_L=$SCREEN_MAX
    if [[ "$N_OBS" -lt "$MAX_L" ]]; then
      MAX_L=$N_OBS
    fi
    generate_group "$GROUP" "$MAX_L" "$SCREEN_STEPS" "$SCREEN_SAMPLES" "screen" "$SCREEN_DIR"
  done
fi

# Also score group 40 matched artifact if present (reference row)
if [[ -f "$OUT_DIR/group_040_generated_matched.nwk" && -f "$OUT_DIR/group_040_generated_matched.fasta" ]]; then
  cp -f "$OUT_DIR/group_040_generated_matched.nwk" "$SCREEN_DIR/group_040_generated_matched.nwk"
  cp -f "$OUT_DIR/group_040_generated_matched.fasta" "$SCREEN_DIR/group_040_generated_matched.fasta"
  cp -f "$OUT_DIR/group_040_observed.nwk" "$SCREEN_DIR/group_040_observed.nwk" 2>/dev/null || true
fi

echo "=== ranking screen gens ==="
$PYTHON scripts/rank_covid_closest_groups.py \
    --data "$DATA" \
    --gen-dir "$SCREEN_DIR" \
    --out-dir "$SCREEN_DIR" \
    --write-md "$SCREEN_DIR/SCREEN_RANKING.md" \
    --promote-top 0

# Pick top PROMOTE_TOP groups excluding 40 from screen ranking
WINNERS=$($PYTHON - <<PY
import json, sys
from pathlib import Path
p = Path("$SCREEN_DIR/closest_groups_ranking.json")
d = json.loads(p.read_text())
ranked = d["best_per_group"]
# Prefer screen-tagged gens for selection; fall back to any non-40
cands = [r for r in ranked if r["group"] != 40 and r.get("tag") == "screen"]
if not cands:
    cands = [r for r in ranked if r["group"] != 40]
top = cands[: int("$PROMOTE_TOP")]
sys.stderr.write(
    "Winners: "
    + str([(r["group"], round(r["combo_distance"], 4), r["n_gen_leaves"]) for r in top])
    + "\n"
)
print(" ".join(str(r["group"]) for r in top))
PY
)
echo "Screen winners: $WINNERS"
if [[ -z "${WINNERS// }" ]]; then
  echo "ERROR: no screen winners"
  exit 1
fi

# ── Phase 2: size-matched regenerate for winners ────────────────────────────
MATCH_DIR="${MATCH_DIR:-results/covid_tree_viz_matched}"
mkdir -p "$MATCH_DIR"

if [[ "$SKIP_MATCH" != "1" ]]; then
  for GROUP in $WINNERS; do
    N_OBS=$(obs_leaves "$GROUP")
    MAX_L=$MATCH_CAP
    if [[ "$N_OBS" -lt "$MAX_L" ]]; then
      MAX_L=$N_OBS
    fi
    generate_group "$GROUP" "$MAX_L" "$MATCH_STEPS" "$MATCH_SAMPLES" "matched" "$MATCH_DIR"
  done
fi

echo "=== ranking matched gens + promoting into $OUT_DIR ==="
$PYTHON scripts/rank_covid_closest_groups.py \
    --data "$DATA" \
    --gen-dir "$MATCH_DIR" \
    --out-dir "$OUT_DIR" \
    --write-md "$OUT_DIR/CLOSEST_GROUPS.md" \
    --promote-top "$PROMOTE_TOP" \
    --exclude-groups 40 \
    --tag-filter matched

# Append group 40 reference row into CLOSEST_GROUPS.md if scorable
if [[ -f "$OUT_DIR/group_040_generated_matched.nwk" ]]; then
  $PYTHON - <<PY
import importlib.util
from pathlib import Path
import sys
sys.path.insert(0, ".")
spec = importlib.util.spec_from_file_location(
    "rank_covid", "scripts/rank_covid_closest_groups.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
from benchmarks.heldout.build_examples import load_tree

out = Path("$OUT_DIR")
obs = load_tree("$DATA", 40)
gen = mod.load_tree_pair(
    out / "group_040_generated_matched.nwk",
    out / "group_040_generated_matched.fasta",
)
m = mod.score_pair(gen, obs, seed=42)
md = out / "CLOSEST_GROUPS.md"
text = md.read_text() if md.exists() else "# Closest COVID gen↔obs groups\n\n"
qd = m["quartet"]
qd_s = f"{qd:.4f}" if qd == qd else "—"
note = (
    "\n## Reference: group 40 (already in viz dir)\n\n"
    "| Group | Leaves | RF | Quartet | Term-edit | Mean best-id | Cov@2% | Combo |\n"
    "|------:|-------:|---:|--------:|----------:|-------------:|-------:|------:|\n"
    f"| 40 | {m['n_gen_leaves']}/{m['n_obs_full']} | {m['rf']:.4f} | "
    f"{qd_s} | {m['terminal_edit']:.4f} | {m['mean_best_identity']:.4f} | "
    f"{m['coverage_eps2pct']:.3f} | {m['combo_distance']:.4f} |\n"
)
if "## Reference: group 40" not in text:
    md.write_text(text.rstrip() + "\n" + note)
print("group40 combo", m["combo_distance"])
PY
fi

# Soft-update README without wiping group 40 notes
$PYTHON - <<PY
from pathlib import Path
import datetime
out = Path("$OUT_DIR")
readme = out / "README.md"
winners = Path("$OUT_DIR/closest_groups_winners.json")
extra = []
if winners.exists():
    import json
    for r in json.loads(winners.read_text()):
        g3 = f"{r['group']:03d}"
        extra.append(
            f"- group {r['group']}: \`group_{g3}_observed.nwk\` ↔ "
            f"\`group_{g3}_generated_matched.nwk\` "
            f"(combo={r['combo_distance']:.4f}, leaves={r['n_gen_leaves']})"
        )
block = (
    "\n\n## Closest additional groups\n\n"
    f"Selected {datetime.datetime.utcnow().strftime('%Y-%m-%d')} from held-out test "
    f"(ckpt={Path('$CKPT').as_posix()}, mrs=$MRS). See \`CLOSEST_GROUPS.md\`.\n\n"
    + ("\n".join(extra) + "\n" if extra else "(none)\n")
)
text = readme.read_text() if readme.exists() else "# COVID tree viz artifacts\n"
if "## Closest additional groups" in text:
    pre = text.split("## Closest additional groups")[0].rstrip()
    text = pre + block
else:
    text = text.rstrip() + block
readme.write_text(text)
print("Updated", readme)
PY

echo "Artifacts:"
ls -la "$OUT_DIR"/group_*_generated_matched.nwk "$OUT_DIR"/CLOSEST_GROUPS.md 2>/dev/null || true
echo "Done: $(date)"
