#!/bin/bash
# Clone VaxSeer + download year-matched lm dominance weights (not full 60GB dump).
#
#   YEAR=2024 SUBTYPE=a_h3n2 bash scripts/setup_vaxseer.sh
#   YEAR=2024 SUBTYPE=a_h1n1 bash scripts/setup_vaxseer.sh

set -euo pipefail
export LABHOME="${LABHOME:-/vast/projects/pranam/lab/nnori}"
ROOT="${VAXSEER_ROOT:-$LABHOME/third_party/vaxseer}"
YEAR="${YEAR:-2024}"
SUBTYPE="${SUBTYPE:-a_h3n2}"
TASK="${TASK:-lm}"

mkdir -p "$(dirname "$ROOT")"
if [[ ! -d "$ROOT/.git" ]]; then
  git clone https://github.com/wxsh1213/vaxseer.git "$ROOT"
else
  echo "VaxSeer already cloned at $ROOT"
fi

cd "$ROOT"
if [[ -f environment.yaml ]]; then
  echo "Optional: conda env create -f environment.yaml"
fi

if [[ -f download_models_from_dropbox.py ]]; then
  python download_models_from_dropbox.py \
    --task "$TASK" \
    --year "$YEAR" \
    --subtype "$SUBTYPE" \
    --output_dir runs
else
  echo "WARN: download_models_from_dropbox.py missing — check upstream README" >&2
  ls -la
fi

echo "Done. VAXSEER_ROOT=$ROOT"
echo "Score with: python scripts/eval_vaxseer_dominance.py --subtype $SUBTYPE --year $YEAR ..."
