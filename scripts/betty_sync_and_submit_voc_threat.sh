#!/bin/bash
# Sync VOC threat Track A artifacts to Betty and submit primary panel jobs.
#
# Usage (laptop, after kinit + Duo-capable SSH):
#   bash scripts/betty_sync_and_submit_voc_threat.sh
#   DRY=1 bash scripts/betty_sync_and_submit_voc_threat.sh   # sync + dry-run submit
#   SUBMIT=0 bash scripts/betty_sync_and_submit_voc_threat.sh  # sync only
#
# Env:
#   BETTY=nnori@login.betty.parcc.upenn.edu
#   REMOTE=~/DiscreteTreeFlows

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

BETTY="${BETTY:-nnori@login.betty.parcc.upenn.edu}"
REMOTE="${REMOTE:-~/DiscreteTreeFlows}"
SUBMIT="${SUBMIT:-1}"
DRY="${DRY:-0}"
CKPT_NAME="${CKPT_NAME:-covid_v7_pmc_hotspot}"
MRS="${MRS:-0.5}"
MAX_LEAVES="${MAX_LEAVES:-250}"
N_STEPS="${N_STEPS:-160}"

echo "=== sync VOC scripts ==="
rsync -avP \
  scripts/voc_threat_lib.py \
  scripts/screen_voc_threat_trees.py \
  scripts/eval_voc_threat_recovery.py \
  scripts/prepare_voc_origin_cases.py \
  scripts/slurm_voc_threat_panel.sh \
  scripts/betty_submit_voc_threat_panel.sh \
  scripts/run_voc_threat_track_a_local.sh \
  scripts/download_sars2_voc_origin_ncbi.py \
  scripts/betty_sync_and_submit_voc_threat.sh \
  scripts/gen_voc_baseline_leaves.py \
  scripts/slurm_voc_threat_baselines.sh \
  "$BETTY:$REMOTE/scripts/"

echo "=== sync panel JSON + results/voc_threat_panel ==="
rsync -avP benchmarks/voc_threat_panel.json "$BETTY:$REMOTE/benchmarks/"
rsync -avP \
  --exclude '*.fasta' --exclude '*.nwk' \
  results/voc_threat_panel/README.md \
  results/voc_threat_panel/cases_manifest_clean.json \
  results/voc_threat_panel/brazil_screen.json \
  results/voc_threat_panel/origin_screen.json \
  results/voc_threat_panel/candidates_test.json \
  results/voc_threat_panel/candidates_train.json \
  results/voc_threat_panel/voc_evescape_annotations.json \
  results/voc_threat_panel/eval_gamma_g003_existing.json \
  "$BETTY:$REMOTE/results/voc_threat_panel/"

# Full primary case folders (observed trees; small)
echo "=== sync primary case folders ==="
python3 - <<'PY'
import json, subprocess, os
from pathlib import Path
root = Path(".")
m = json.loads((root / "results/voc_threat_panel/cases_manifest_clean.json").read_text())
betty = os.environ.get("BETTY", "nnori@login.betty.parcc.upenn.edu")
remote = os.environ.get("REMOTE", "~/DiscreteTreeFlows")
for c in m["cases"]:
    local = root / c["case_dir"]
    if not local.exists():
        print("SKIP missing", local)
        continue
    dest = f"{betty}:{remote}/{c['case_dir']}/"
    subprocess.run(["ssh", betty, f"mkdir -p {remote}/{c['case_dir']}"], check=True)
    subprocess.run(
        ["rsync", "-avP",
         str(local / "case.json"),
         str(local / "observed_anc_aa.fasta"),
         str(local / "observed.nwk"),
         dest],
        check=False,
    )
    meta = local / "meta.csv"
    if meta.exists():
        subprocess.run(["rsync", "-avP", str(meta), dest], check=False)
print("synced", len(m["cases"]), "cases")
PY

if [[ "$SUBMIT" != "1" ]]; then
  echo "Sync only (SUBMIT=0)"
  exit 0
fi

echo "=== submit on Betty (DRY=$DRY) ==="
ssh "$BETTY" bash -s <<REMOTE_EOF
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:\$PATH
# Slurm is only on PATH for login shells; \`ssh betty bash -s\` is not one.
export PATH=/cm/local/apps/slurm/current/bin:/vast/parcc/sw/bin:\$PATH
export LABHOME=\${LABHOME:-/vast/projects/pranam/lab/nnori}
cd ~/DiscreteTreeFlows
mkdir -p logs results/voc_threat_panel

# Ensure ckpt exists (lab backup fallback)
CKPT=checkpoints/${CKPT_NAME}/best.pt
if [[ ! -f "\$CKPT" ]]; then
  B=\$LABHOME/checkpoints_backup/${CKPT_NAME}/best.pt
  if [[ -f "\$B" ]]; then
    mkdir -p checkpoints/${CKPT_NAME}
    cp -n "\$B" "\$CKPT" || cp "\$B" "\$CKPT"
  fi
fi
ls -la "\$CKPT" || { echo "ERROR: missing \$CKPT"; exit 1; }

# Sanity: primary groups exist on Betty
python3 - <<'PY'
import json
from pathlib import Path
m=json.loads(Path("results/voc_threat_panel/cases_manifest_clean.json").read_text())
missing=[]
for c in m["cases"]:
    split, gid = c["split"], int(c["gid"])
    p=Path(f"data/covid/{split}/group_{gid:03d}_anc_aa.fasta")
    if not p.exists():
        missing.append(str(p))
if missing:
    raise SystemExit("Missing on Betty:\\n" + "\\n".join(missing))
print(f"OK {len(m['cases'])} primary cases have anc_aa on Betty")
PY

DRY=${DRY} CKPT_NAME=${CKPT_NAME} MRS=${MRS} MAX_LEAVES=${MAX_LEAVES} N_STEPS=${N_STEPS} \\
  bash scripts/betty_submit_voc_threat_panel.sh

# Also submit NeutralBD baselines for each primary case (parallel, no dependency)
echo "=== submit NeutralBD baselines ==="
python3 - <<'PY'
import json, os, subprocess
from pathlib import Path
m = json.loads(Path("results/voc_threat_panel/cases_manifest_clean.json").read_text())
dry = os.environ.get("DRY", "0") == "1"
jobs = []
for c in m["cases"]:
    env = {
        **os.environ,
        "CASE_DIR": c["case_dir"],
        "VOC_ID": c["voc_id"],
        "N_LEAVES": "250",
        "METHODS": "neutral_bd",
        "DEVICE": "cuda",
    }
    cmd = ["sbatch", "--export=ALL,CASE_DIR,VOC_ID,N_LEAVES,METHODS,DEVICE",
           "scripts/slurm_voc_threat_baselines.sh"]
    print(("[DRY] " if dry else "") + f"baseline {c['voc_id']} {c['case_dir']}")
    if dry:
        jobs.append({"voc_id": c["voc_id"], "job_id": None})
        continue
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, check=True)
    jid = next((int(t) for t in r.stdout.split() if t.isdigit()), None)
    print(r.stdout.strip())
    jobs.append({"voc_id": c["voc_id"], "case_dir": c["case_dir"], "job_id": jid})
Path("results/voc_threat_panel/baseline_submit_ids.json").write_text(
    json.dumps({"jobs": jobs}, indent=2) + "\n"
)
print("wrote results/voc_threat_panel/baseline_submit_ids.json")
PY

echo "=== squeue ==="
squeue -u nnori | head -40
REMOTE_EOF

echo "=== done local side ==="
