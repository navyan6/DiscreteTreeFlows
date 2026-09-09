#!/bin/bash
# Submit VOC threat-panel generation+eval jobs on Betty from cases_manifest_clean.json.
#
# Prerequisites (local or Betty):
#   python scripts/screen_voc_threat_trees.py --data-dir data/covid/test --split test ...
#   python scripts/screen_voc_threat_trees.py --data-dir data/covid/train --split train --groups ...
#   python scripts/prepare_voc_origin_cases.py --candidates ... --case Gamma:test:3 ...
#
# Usage (on Betty):
#   bash scripts/betty_submit_voc_threat_panel.sh
#   DRY=1 bash scripts/betty_submit_voc_threat_panel.sh

set -euo pipefail
cd "${HOME}/DiscreteTreeFlows" 2>/dev/null || cd "$(dirname "$0")/.."

# Betty's login nodes only put Slurm on PATH for login shells, so a plain
# `ssh betty bash script.sh` cannot find sbatch. Add it explicitly.
export PATH=/cm/local/apps/slurm/current/bin:/vast/parcc/sw/bin:$PATH

MANIFEST="${MANIFEST:-results/voc_threat_panel/cases_manifest_clean.json}"
CKPT_NAME="${CKPT_NAME:-covid_v7_pmc_hotspot}"
MRS="${MRS:-0.5}"
MAX_LEAVES="${MAX_LEAVES:-250}"
N_STEPS="${N_STEPS:-160}"
BRANCH_SCALE="${BRANCH_SCALE:-6.0}"
SEED="${SEED:-42}"
DRY="${DRY:-0}"
OUT_IDS="${OUT_IDS:-results/voc_threat_panel/submit_ids.json}"

if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: missing $MANIFEST — run scripts/run_voc_threat_track_a_local.sh first"
  exit 1
fi
echo "Using MANIFEST=$MANIFEST"

mkdir -p logs results/voc_threat_panel

export MANIFEST CKPT_NAME MRS MAX_LEAVES N_STEPS BRANCH_SCALE SEED DRY OUT_IDS

python3 - <<'PY'
import json, os, subprocess
from pathlib import Path

manifest = json.loads(Path(os.environ["MANIFEST"]).read_text())
cases = manifest["cases"]
dry = os.environ.get("DRY", "0") == "1"
ckpt = os.environ["CKPT_NAME"]
mrs = os.environ["MRS"]
max_leaves = os.environ["MAX_LEAVES"]
n_steps = os.environ["N_STEPS"]
branch = os.environ["BRANCH_SCALE"]
seed = os.environ["SEED"]
out_ids = Path(os.environ["OUT_IDS"])

jobs = []
for c in cases:
    case_dir = c["case_dir"]
    voc = c["voc_id"]
    split = c["split"]
    gid = int(c["gid"])
    data = c.get("source_data_dir") or f"data/covid/{split}"
    # normalize to repo-relative
    if data.startswith("/"):
        # keep basename path under data/covid
        data = f"data/covid/{split}"

    env = {
        **os.environ,
        "CASE_DIR": case_dir,
        "VOC_ID": voc,
        "DATA": data,
        "GROUP": str(gid),
        "CKPT_NAME": ckpt,
        "MRS": mrs,
        "MAX_LEAVES": max_leaves,
        "N_STEPS": n_steps,
        "BRANCH_SCALE": branch,
        "SEED": seed,
        "REF_TIP": "",
    }
    if voc == "Gamma" and gid == 3:
        env["REF_TIP"] = "MZ397166.1"

    export_list = (
        "ALL,CASE_DIR,VOC_ID,DATA,GROUP,CKPT_NAME,MRS,MAX_LEAVES,"
        "N_STEPS,BRANCH_SCALE,SEED,REF_TIP"
    )
    cmd = ["sbatch", f"--export={export_list}", "scripts/slurm_voc_threat_panel.sh"]
    print(f"{'[DRY] ' if dry else ''}submit {voc} {split} g{gid:03d} -> {case_dir}")
    if dry:
        jobs.append({"voc_id": voc, "case_dir": case_dir, "data": data, "group": gid, "job_id": None})
        continue
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, check=True)
    jid = None
    for tok in r.stdout.split():
        if tok.isdigit():
            jid = int(tok)
    print(r.stdout.strip())
    jobs.append({
        "voc_id": voc,
        "case_dir": case_dir,
        "data": data,
        "group": gid,
        "job_id": jid,
        "stdout": r.stdout.strip(),
    })

out_ids.parent.mkdir(parents=True, exist_ok=True)
out_ids.write_text(json.dumps({
    "ckpt": ckpt,
    "mrs": float(mrs),
    "max_leaves": int(max_leaves),
    "n_steps": int(n_steps),
    "jobs": jobs,
}, indent=2) + "\n")
print(f"wrote {out_ids} ({len(jobs)} jobs)")
PY
