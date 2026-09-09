#!/bin/bash
#SBATCH --job-name=h3n2_treeviz
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/h3n2_tree_viz_%j.log
#SBATCH --error=logs/h3n2_tree_viz_%j.log
#
# Screen a few held-out H3N2 test groups, rank gen↔obs closeness, then
# regenerate size-matched trees for the top-2 winners into results/h3n2_tree_viz/.
#
# Env overrides:
#   CAND_GROUPS="5 20 40 46 55 59"
#   SCREEN_MAX=96 SCREEN_STEPS=60 SCREEN_SAMPLES=1
#   MATCH_CAP=250 MATCH_STEPS=160 MATCH_SAMPLES=1 PROMOTE_TOP=2
#   MRS=0.5 CKPT_NAME=h3n2_v3_lit_hotspot
#
# Usage:
#   sbatch scripts/slurm_h3n2_tree_viz.sh

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
export TRANSFORMERS_CACHE=${TRANSFORMERS_CACHE:-$HF_HOME}
export TORCH_HOME=${TORCH_HOME:-$LABHOME/torch_cache}
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs results/h3n2_tree_viz results/h3n2_tree_viz_screen results/h3n2_tree_viz_matched

CKPT_NAME="${CKPT_NAME:-h3n2_v3_lit_hotspot}"
CKPT="checkpoints/${CKPT_NAME}/best.pt"
DATA="${DATA:-data/h3n2/test}"
MRS="${MRS:-0.5}"
BRANCH_SCALE="${BRANCH_SCALE:-6.0}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"

# Temporally diverse held-out groups (all ~378–400 tips). Smallest first.
CAND_GROUPS="${CAND_GROUPS:-46 5 20 40 55 59}"
SCREEN_MAX="${SCREEN_MAX:-96}"
SCREEN_STEPS="${SCREEN_STEPS:-60}"
SCREEN_SAMPLES="${SCREEN_SAMPLES:-1}"
MATCH_CAP="${MATCH_CAP:-250}"
MATCH_STEPS="${MATCH_STEPS:-160}"
MATCH_SAMPLES="${MATCH_SAMPLES:-1}"
PROMOTE_TOP="${PROMOTE_TOP:-2}"
SKIP_SCREEN="${SKIP_SCREEN:-0}"
SKIP_MATCH="${SKIP_MATCH:-0}"

SCREEN_DIR="${SCREEN_DIR:-results/h3n2_tree_viz_screen}"
MATCH_DIR="${MATCH_DIR:-results/h3n2_tree_viz_matched}"
OUT_DIR="${OUT_DIR:-results/h3n2_tree_viz}"
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

group_dates() {
  local g="$1"
  $PYTHON - <<PY
import csv
from pathlib import Path
p = Path("$DATA") / f"h3n2test_group_{int('$g'):03d}.csv"
if not p.exists():
    print("?")
else:
    dates = [r["date"] for r in csv.DictReader(p.open()) if r.get("date")]
    print(f"{min(dates)}..{max(dates)}" if dates else "?")
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

EPOCH=$($PYTHON - <<PY
import torch
c=torch.load("$CKPT", map_location="cpu", weights_only=False)
print(c.get("epoch", "?"))
PY
)

echo "Start: $(date)"
echo "ckpt=$CKPT epoch=$EPOCH mrs=$MRS cand_groups=[$CAND_GROUPS]"
echo "screen: max=$SCREEN_MAX steps=$SCREEN_STEPS samples=$SCREEN_SAMPLES"
echo "match:  cap=$MATCH_CAP steps=$MATCH_STEPS samples=$MATCH_SAMPLES promote=$PROMOTE_TOP"

# ── Phase 1: screen ─────────────────────────────────────────────────────────
if [[ "$SKIP_SCREEN" != "1" ]]; then
  for GROUP in $CAND_GROUPS; do
    if [[ ! -f "$DATA/group_$(printf '%03d' "$GROUP")_rooted.nwk" ]]; then
      echo "WARN: skip missing group $GROUP"
      continue
    fi
    N_OBS=$(obs_leaves "$GROUP")
    MAX_L=$SCREEN_MAX
    if [[ "$N_OBS" -lt "$MAX_L" ]]; then
      MAX_L=$N_OBS
    fi
    generate_group "$GROUP" "$MAX_L" "$SCREEN_STEPS" "$SCREEN_SAMPLES" "screen" "$SCREEN_DIR"
  done
fi

echo "=== ranking screen gens ==="
$PYTHON scripts/rank_covid_closest_groups.py \
    --data "$DATA" \
    --gen-dir "$SCREEN_DIR" \
    --out-dir "$SCREEN_DIR" \
    --write-md "$SCREEN_DIR/SCREEN_RANKING.md" \
    --promote-top 0

WINNERS=$($PYTHON - <<PY
import json, sys
from pathlib import Path
p = Path("$SCREEN_DIR/closest_groups_ranking.json")
d = json.loads(p.read_text())
ranked = d["best_per_group"]
cands = [r for r in ranked if r.get("tag") == "screen"]
if not cands:
    cands = list(ranked)
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
    --tag-filter matched

# Ensure observed_matched exists for every winner (size-matched tip prune)
$PYTHON - <<PY
import importlib.util
import json
from pathlib import Path
import sys
sys.path.insert(0, ".")
from benchmarks.heldout.build_examples import load_tree
from benchmarks.metrics import trees as T

spec = importlib.util.spec_from_file_location(
    "rank_covid", "scripts/rank_covid_closest_groups.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

out = Path("$OUT_DIR")
winners = json.loads((out / "closest_groups_winners.json").read_text())
for r in winners:
    gid = r["group"]
    g3 = f"{gid:03d}"
    obs = load_tree("$DATA", gid)
    gen = mod.load_tree_pair(Path(r["gen_nwk"]), Path(r["gen_fasta"]))
    n_gen = len(T.leaf_labels(gen))
    obs_m = mod.size_match_obs(obs, n_gen, seed=42)
    (out / f"group_{g3}_observed_matched.nwk").write_text(T.to_newick(obs_m) + "\n")
    print(f"observed_matched group={gid} leaves={n_gen}")
PY

# README (quoted heredoc — backticks/fences must not be shell-expanded)
export H3N2_VIZ_OUT_DIR="$OUT_DIR" H3N2_VIZ_DATA="$DATA" H3N2_VIZ_CKPT="$CKPT" \
  H3N2_VIZ_CKPT_NAME="$CKPT_NAME" H3N2_VIZ_EPOCH="$EPOCH" H3N2_VIZ_MRS="$MRS" \
  H3N2_VIZ_MATCH_CAP="$MATCH_CAP" H3N2_VIZ_MATCH_STEPS="$MATCH_STEPS" \
  H3N2_VIZ_SCREEN_MAX="$SCREEN_MAX" H3N2_VIZ_SCREEN_STEPS="$SCREEN_STEPS" \
  H3N2_VIZ_BRANCH_SCALE="$BRANCH_SCALE" H3N2_VIZ_MAX_SEQ_LEN="$MAX_SEQ_LEN"
$PYTHON - <<'PY'
import csv, datetime, json, os
from pathlib import Path
from ete3 import Tree

out = Path(os.environ["H3N2_VIZ_OUT_DIR"])
data = Path(os.environ["H3N2_VIZ_DATA"])
winners = json.loads((out / "closest_groups_winners.json").read_text())
rows = []
for r in winners:
    gid = r["group"]
    g3 = f"{gid:03d}"
    csvp = data / f"h3n2test_group_{g3}.csv"
    d0 = d1 = "?"
    if csvp.exists():
        dates = [x["date"] for x in csv.DictReader(csvp.open()) if x.get("date")]
        if dates:
            d0, d1 = min(dates), max(dates)
    obs_n = len(Tree(str(out / f"group_{g3}_observed.nwk"), format=1).get_leaves())
    om = out / f"group_{g3}_observed_matched.nwk"
    gm = out / f"group_{g3}_generated_matched.nwk"
    om_n = len(Tree(str(om), format=1).get_leaves()) if om.exists() else "?"
    gm_n = len(Tree(str(gm), format=1).get_leaves()) if gm.exists() else "?"
    rows.append((gid, g3, d0, d1, obs_n, om_n, gm_n, r["combo_distance"]))

job = os.environ.get("SLURM_JOB_ID", "local")
ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%MZ")
lines = [
    "# H3N2 tree viz artifacts",
    "",
    "Paired **observed_matched** vs **generated_matched** Newick trees for held-out H3N2 test groups.",
    "Same root/group; observed tip-subsampled to generated leaf count for side-by-side viz.",
    "",
    "| Field | Value |",
    "|-------|-------|",
    f"| Split | `{os.environ['H3N2_VIZ_DATA']}` (held-out test) |",
    f"| Checkpoint | `{os.environ['H3N2_VIZ_CKPT']}` (`{os.environ['H3N2_VIZ_CKPT_NAME']}`, epoch {os.environ['H3N2_VIZ_EPOCH']}) |",
    f"| Mutation rate scale (mrs) | {os.environ['H3N2_VIZ_MRS']} |",
    f"| Match max_leaves / n_steps | {os.environ['H3N2_VIZ_MATCH_CAP']} / {os.environ['H3N2_VIZ_MATCH_STEPS']} |",
    f"| Screen max_leaves / n_steps | {os.environ['H3N2_VIZ_SCREEN_MAX']} / {os.environ['H3N2_VIZ_SCREEN_STEPS']} |",
    f"| branch_rate_scale | {os.environ['H3N2_VIZ_BRANCH_SCALE']} |",
    f"| max_seq_len | {os.environ['H3N2_VIZ_MAX_SEQ_LEN']} |",
    f"| Job | SLURM {job} @ {ts} |",
    "",
    "## Selected groups (closest gen↔obs among screen candidates)",
    "",
    "| Group | Dates | Full obs leaves | Matched leaves | Combo ↓ | Viz pair |",
    "|------:|-------|----------------:|---------------:|--------:|----------|",
]
for gid, g3, d0, d1, obs_n, om_n, gm_n, combo in rows:
    lines.append(
        f"| {gid} | {d0}..{d1} | {obs_n} | {om_n}/{gm_n} | {combo:.4f} | "
        f"`group_{g3}_observed_matched.nwk` ↔ `group_{g3}_generated_matched.nwk` |"
    )
lines += [
    "",
    "Location: not in local CSV (name/date only); sequences are GISAID H3N2 HA.",
    "",
    "## Files",
    "",
    "- `group_XXX_observed.nwk` — full ground-truth rooted tree",
    "- `group_XXX_observed_matched.nwk` — observed pruned to gen leaf count (seed=42)",
    "- `group_XXX_generated_matched.nwk` — size-matched model sample",
    "- `CLOSEST_GROUPS.md` / ranking CSV+JSON — closeness scores",
    "",
]
(out / "README.md").write_text("\n".join(lines) + "\n")
print("Wrote", out / "README.md")
for gid, g3, *_ in rows:
    print("PAIR", out / f"group_{g3}_observed_matched.nwk")
    print("PAIR", out / f"group_{g3}_generated_matched.nwk")
PY

# Sync local README if present is fine; Betty copy now exists too
echo "Artifacts in $OUT_DIR:"
ls -la "$OUT_DIR"
echo "Done: $(date)"
