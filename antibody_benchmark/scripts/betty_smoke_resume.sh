#!/usr/bin/env bash
# Resume incremental model smokes on Betty after SSH/Kerberos is back.
set -euo pipefail
cd "${HOME}/DiscreteTreeFlows"
PY=/vast/projects/pranam/lab/nnori/.conda/envs/treesbm/bin/python
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
export AB_DASM_DIR="${AB_DASM_DIR:-$HOME/antibody_benchmark_raw/dasm/extracted/dasm-experiments-data}"

# Fix incomplete thrifty/dasm pretrained extracts (cached zip but missing .yml)
NETAM_PRE="$PWD/antibody_benchmark/data/raw/repos/netam/netam/_pretrained"
for zip in thrifty-0.2.0.zip dasm-1.0.0.zip; do
  if [[ -f "$NETAM_PRE/$zip" ]]; then
    # Re-extract if expected model yml missing
    need=0
    if [[ "$zip" == thrifty-0.2.0.zip && ! -f "$NETAM_PRE/thrifty-models-0.2.0/models/ThriftyHumV0.2-59.yml" ]]; then need=1; fi
    if [[ "$zip" == dasm-1.0.0.zip && ! -f "$NETAM_PRE/dasm-models-1.0.0/models/DASMHumV1.0-4M.yml" ]]; then need=1; fi
    if [[ $need -eq 1 ]]; then
      echo "Re-extracting $zip"
      rm -rf "$NETAM_PRE/${zip%.zip}"* 2>/dev/null || true
      # keep zip; remove extracted dir only
      if [[ "$zip" == thrifty-0.2.0.zip ]]; then rm -rf "$NETAM_PRE/thrifty-models-0.2.0"; fi
      if [[ "$zip" == dasm-1.0.0.zip ]]; then rm -rf "$NETAM_PRE/dasm-models-1.0.0"; fi
      unzip -o "$NETAM_PRE/$zip" -d "$NETAM_PRE"
    fi
  fi
done

mkdir -p antibody_benchmark/data/raw/cosine/checkpoints antibody_benchmark/results/smoke_reports

# 1) Thrifty
$PY antibody_benchmark/scripts/smoke_models.py --models thrifty --family-index 0
# 2) DASM+Thrifty
$PY antibody_benchmark/scripts/smoke_models.py --models dasm_thrifty --family-index 0
# 3) CoSiNE (needs full ckpt ~1.7G)
CKPT=antibody_benchmark/data/raw/cosine/checkpoints/cosine_dasm.ckpt
if [[ -f "$CKPT" ]]; then
  sz=$(stat -c%s "$CKPT" 2>/dev/null || stat -f%z "$CKPT")
  if [[ "$sz" -lt 1000000000 ]]; then
    echo "BLOCKER: cosine ckpt incomplete ($sz bytes); scp full 1.7G file first"
  else
    # ensure cosine importable
    pip_bin=/vast/projects/pranam/lab/nnori/.conda/envs/treesbm/bin/pip
    $pip_bin install -e antibody_benchmark/data/raw/repos/cosine 2>&1 | tail -20 || true
    $PY antibody_benchmark/scripts/smoke_models.py --models cosine --family-index 0
  fi
else
  echo "BLOCKER: missing $CKPT"
fi
# 4) TreeSBM
$PY antibody_benchmark/scripts/smoke_models.py --models treesbm --family-index 0

echo '--- scale: 1 fam x 1 already done above; next 5 fam x 2 ---'
# Optional scale after one-tree passes — edit enabled models in a temp config
