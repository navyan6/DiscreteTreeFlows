#!/bin/bash
#SBATCH --job-name=h1n1_treeviz
#SBATCH --partition=b200-mig45
#SBATCH --qos=mig-max
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/h1n1_tree_viz_%j.log
#SBATCH --error=logs/h1n1_tree_viz_%j.log
#
# Size-matched obs vs gen Newicks for 2 held-out H1N1 groups from the SAME
# region/population (from data/h1n1_temporal/SPLIT_PROTOCOL.json units).
#
# Default: latest temporal ckpt + USA: Michigan groups 9 and 13.
#
# Env overrides:
#   CKPT_NAME=h1n1_temporal_v1
#   DATA=data/h1n1_temporal/test
#   TREE_GROUPS="9 13"   (do NOT use GROUPS — Slurm/OS sets that)
#   REGION="USA: Michigan"
#   MRS=0.5 MATCH_CAP=250 MATCH_STEPS=160
#
# Usage:
#   sbatch scripts/slurm_h1n1_tree_viz.sh

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
export LABHOME=${LABHOME:-/vast/projects/pranam/lab/nnori}
export HF_HOME=${HF_HOME:-$LABHOME/hf_cache}
export TRANSFORMERS_CACHE=${TRANSFORMERS_CACHE:-$HF_HOME}
export TORCH_HOME=${TORCH_HOME:-$LABHOME/torch_cache}
export PYTHONUNBUFFERED=1

cd ~/DiscreteTreeFlows
mkdir -p logs results/h1n1_tree_viz

CKPT_NAME="${CKPT_NAME:-h1n1_temporal_v1}"
CKPT="checkpoints/${CKPT_NAME}/best.pt"
DATA="${DATA:-data/h1n1_temporal/test}"
MRS="${MRS:-0.5}"
BRANCH_SCALE="${BRANCH_SCALE:-6.0}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-566}"
# Avoid env name GROUPS (often pre-set by the OS / Slurm to numeric GIDs).
TREE_GROUPS="${TREE_GROUPS:-9 13}"
REGION="${REGION:-USA: Michigan}"
MATCH_CAP="${MATCH_CAP:-250}"
MATCH_STEPS="${MATCH_STEPS:-160}"
MATCH_SAMPLES="${MATCH_SAMPLES:-1}"
OUT_DIR="${OUT_DIR:-results/h1n1_tree_viz}"
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
echo "ckpt=$CKPT epoch=$EPOCH mrs=$MRS region=$REGION groups=[$TREE_GROUPS]"
echo "match: cap=$MATCH_CAP steps=$MATCH_STEPS samples=$MATCH_SAMPLES"

for GROUP in $TREE_GROUPS; do
  G3=$(printf '%03d' "$GROUP")
  if [[ ! -f "$DATA/group_${G3}_rooted.nwk" ]]; then
    echo "ERROR: missing $DATA/group_${G3}_rooted.nwk"
    exit 1
  fi
  N_OBS=$(obs_leaves "$GROUP")
  MAX_L=$MATCH_CAP
  if [[ "$N_OBS" -lt "$MAX_L" ]]; then
    MAX_L=$N_OBS
  fi
  generate_group "$GROUP" "$MAX_L" "$MATCH_STEPS" "$MATCH_SAMPLES" "matched" "$OUT_DIR"
done

# Size-match observed tips to generated leaf count
$PYTHON - <<PY
import importlib.util
import sys
from pathlib import Path
sys.path.insert(0, ".")
from benchmarks.heldout.build_examples import load_tree
from benchmarks.metrics import trees as T

spec = importlib.util.spec_from_file_location(
    "rank_covid", "scripts/rank_covid_closest_groups.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

out = Path("$OUT_DIR")
groups = [int(x) for x in "$TREE_GROUPS".split()]
for gid in groups:
    g3 = f"{gid:03d}"
    obs = load_tree("$DATA", gid)
    gen_nwk = out / f"group_{g3}_generated_matched.nwk"
    gen_fa = out / f"group_{g3}_generated_matched.fasta"
    gen = mod.load_tree_pair(gen_nwk, gen_fa if gen_fa.exists() else None)
    n_gen = len(T.leaf_labels(gen))
    obs_m = mod.size_match_obs(obs, n_gen, seed=42)
    (out / f"group_{g3}_observed_matched.nwk").write_text(T.to_newick(obs_m) + "\n")
    # also keep a primary alias without seed suffix
    print(f"observed_matched group={gid} leaves={n_gen}")
PY

export H1N1_VIZ_OUT_DIR="$OUT_DIR" H1N1_VIZ_DATA="$DATA" H1N1_VIZ_CKPT="$CKPT" \
  H1N1_VIZ_CKPT_NAME="$CKPT_NAME" H1N1_VIZ_EPOCH="$EPOCH" H1N1_VIZ_MRS="$MRS" \
  H1N1_VIZ_MATCH_CAP="$MATCH_CAP" H1N1_VIZ_MATCH_STEPS="$MATCH_STEPS" \
  H1N1_VIZ_BRANCH_SCALE="$BRANCH_SCALE" H1N1_VIZ_MAX_SEQ_LEN="$MAX_SEQ_LEN" \
  H1N1_VIZ_GROUPS="$TREE_GROUPS" H1N1_VIZ_REGION="$REGION"
$PYTHON - <<'PY'
import csv, datetime, json, os
from pathlib import Path
from ete3 import Tree

out = Path(os.environ["H1N1_VIZ_OUT_DIR"])
data = Path(os.environ["H1N1_VIZ_DATA"])
region = os.environ["H1N1_VIZ_REGION"]
groups = [int(x) for x in os.environ["H1N1_VIZ_GROUPS"].split()]

# Prefer SPLIT_PROTOCOL unit mapping when present
proto_path = data.parent / "SPLIT_PROTOCOL.json"
unit_by_g = {}
if proto_path.exists():
    proto = json.loads(proto_path.read_text())
    for u in proto.get("units", []):
        for g in u.get("groups", []):
            unit_by_g[int(g)] = u.get("unit", "?")

rows = []
for gid in groups:
    g3 = f"{gid:03d}"
    meta = data / f"group_{g3}_meta.csv"
    d0 = d1 = "?"
    if meta.exists():
        dates = [x["date"] for x in csv.DictReader(meta.open()) if x.get("date")]
        if dates:
            d0, d1 = min(dates), max(dates)
    obs_n = len(Tree(str(out / f"group_{g3}_observed.nwk"), format=1).get_leaves())
    om = out / f"group_{g3}_observed_matched.nwk"
    gm = out / f"group_{g3}_generated_matched.nwk"
    om_n = len(Tree(str(om), format=1).get_leaves()) if om.exists() else "?"
    gm_n = len(Tree(str(gm), format=1).get_leaves()) if gm.exists() else "?"
    unit = unit_by_g.get(gid, region)
    rows.append((gid, g3, unit, d0, d1, obs_n, om_n, gm_n))

job = os.environ.get("SLURM_JOB_ID", "local")
ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%MZ")
lines = [
    "# H1N1 tree viz artifacts",
    "",
    "Paired **observed_matched** vs **generated_matched** Newick trees for 2 held-out",
    "H1N1 test groups from the **same region/population**.",
    "",
    "| Field | Value |",
    "|-------|-------|",
    f"| Region / unit | `{region}` |",
    f"| Split | `{os.environ['H1N1_VIZ_DATA']}` (temporal held-out test, year 2025) |",
    f"| Checkpoint | `{os.environ['H1N1_VIZ_CKPT']}` (`{os.environ['H1N1_VIZ_CKPT_NAME']}`, epoch {os.environ['H1N1_VIZ_EPOCH']}) |",
    f"| Mutation rate scale (mrs) | {os.environ['H1N1_VIZ_MRS']} |",
    f"| Match max_leaves / n_steps | {os.environ['H1N1_VIZ_MATCH_CAP']} / {os.environ['H1N1_VIZ_MATCH_STEPS']} |",
    f"| branch_rate_scale | {os.environ['H1N1_VIZ_BRANCH_SCALE']} |",
    f"| max_seq_len | {os.environ['H1N1_VIZ_MAX_SEQ_LEN']} |",
    f"| Job | SLURM {job} @ {ts} |",
    "",
    "## Selected groups (same region)",
    "",
    "| Group | Unit | Dates | Full obs leaves | Matched leaves | Viz pair |",
    "|------:|------|-------|----------------:|---------------:|----------|",
]
for gid, g3, unit, d0, d1, obs_n, om_n, gm_n in rows:
    lines.append(
        f"| {gid} | {unit} | {d0}..{d1} | {obs_n} | {om_n}/{gm_n} | "
        f"`group_{g3}_observed_matched.nwk` ↔ `group_{g3}_generated_matched.nwk` |"
    )
lines += [
    "",
    "## Notes",
    "",
    "- Latest H1N1 run is **temporal** (`h1n1_temporal_v1`), not geo-holdout.",
    "- Groups are still built per location/country unit (see `SPLIT_PROTOCOL.json`),",
    "  so same-region pairs are available within the temporal test split.",
    "- Per-group `meta.csv` has name/date only; region comes from the protocol unit map.",
    "",
    "## Recommended viz pairs (equal leaf count)",
    "",
]
for gid, g3, *_ in rows:
    lines.append(
        f"{len([x for x in rows if x[0] <= gid])}. **group {gid}:** "
        f"`group_{g3}_observed_matched.nwk` ↔ `group_{g3}_generated_matched.nwk`"
    )
lines += [
    "",
    "## Quick viz",
    "",
    "```bash",
    "python - <<'PY'",
    "from ete3 import Tree",
    f"for g in {tuple(f'{gid:03d}' for gid in groups)}:",
    '    obs = Tree(f"group_{g}_observed_matched.nwk", format=1)',
    '    gen = Tree(f"group_{g}_generated_matched.nwk", format=1)',
    '    print(g, "leaves", len(obs), len(gen))',
    "PY",
    "```",
    "",
    "## Regenerate on Betty",
    "",
    "```bash",
    f'TREE_GROUPS="{" ".join(str(g) for g in groups)}" REGION="{region}" \\',
    "  sbatch scripts/slurm_h1n1_tree_viz.sh",
    "```",
    "",
]
(out / "README.md").write_text("\n".join(lines) + "\n")
print("Wrote", out / "README.md")
for gid, g3, *_ in rows:
    print("PAIR", out / f"group_{g3}_observed_matched.nwk")
    print("PAIR", out / f"group_{g3}_generated_matched.nwk")
PY

echo "Artifacts in $OUT_DIR:"
ls -la "$OUT_DIR"
echo "Done: $(date)"
