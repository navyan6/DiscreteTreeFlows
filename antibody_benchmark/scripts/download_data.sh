#!/usr/bin/env bash
# Download official antibody benchmark assets.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RAW="$ROOT/data/raw"
mkdir -p "$RAW"/{dasm,cosine/checkpoints,thrifty,repos}

echo "== CoSiNE HF dataset =="
hf download --repo-type dataset thematrixmaster/cosine --local-dir "$RAW/cosine/hf_dataset"

echo "== CoSiNE checkpoint =="
hf download thematrixmaster/cosine cosine_dasm.ckpt --local-dir "$RAW/cosine/checkpoints"

echo "== DASM Zenodo (may 403 on some networks; try Betty) =="
ZURL="https://zenodo.org/api/records/17322891/files/dasm-experiments-data.tar.gz/content"
if curl -L --fail -o "$RAW/dasm/dasm-experiments-data.tar.gz" "$ZURL"; then
  mkdir -p "$RAW/dasm/extracted"
  tar -xzf "$RAW/dasm/dasm-experiments-data.tar.gz" -C "$RAW/dasm/extracted"
  echo "DASM extracted to $RAW/dasm/extracted"
else
  echo "BLOCKER: Zenodo download failed. On Betty:"
  echo "  curl -L --fail -o ~/antibody_benchmark_raw/dasm/dasm-experiments-data.tar.gz '$ZURL'"
  echo "  scp betty:~/antibody_benchmark_raw/dasm/dasm-experiments-data.tar.gz $RAW/dasm/"
fi

echo "== Thrifty Dryad (often 401 without browser auth) =="
DRYAD="https://datadryad.org/api/v2/files/4098963/download"
if curl -L --fail -o "$RAW/thrifty/for_dryad.v2.zip" "$DRYAD"; then
  mkdir -p "$RAW/thrifty/extracted"
  unzip -o "$RAW/thrifty/for_dryad.v2.zip" -d "$RAW/thrifty/extracted"
else
  echo "BLOCKER: Dryad API returned auth error. Download manually from"
  echo "  https://doi.org/10.5061/dryad.np5hqc044  (for_dryad.v2.zip)"
fi

echo "== Clone official repos (optional) =="
cd "$RAW/repos"
for repo in songlab-cal/cosine matsengrp/netam matsengrp/dasm-experiments matsengrp/thrifty-experiments-1; do
  name="${repo##*/}"
  if [[ ! -d "$name" ]]; then
    git clone --depth 1 "https://github.com/${repo}.git" "$name" || true
  fi
done

echo "Done. See antibody_benchmark/DATA_AUDIT.md"
