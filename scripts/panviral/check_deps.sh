#!/usr/bin/env bash
# Fail fast if this machine cannot run the pan-viral data pipeline.
#
#   bash scripts/panviral/check_deps.sh
#   TREESBM_PY=/path/to/python bash scripts/panviral/check_deps.sh
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"

resolve_py() {
    if [ -n "${TREESBM_PY:-}" ]; then
        echo "$TREESBM_PY"
        return
    fi
    if [ -x "${HOME}/.conda/envs/treesbm-data/bin/python" ]; then
        echo "${HOME}/.conda/envs/treesbm-data/bin/python"
        return
    fi
    if [ -x "${HOME}/.conda/envs/treesbm/bin/python" ]; then
        echo "${HOME}/.conda/envs/treesbm/bin/python"
        return
    fi
    command -v python3 || command -v python
}

PY="$(resolve_py)"
BIN="$(dirname "$PY")"
export PATH="$BIN:$PATH"

echo "python: $PY"
"$PY" - <<'PY'
import sys
try:
    import Bio
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError as e:
    print("MISSING python package:", e, file=sys.stderr)
    print("  conda env create -f scripts/panviral/environment.yml", file=sys.stderr)
    print("  # or: pip install -r scripts/panviral/requirements.txt", file=sys.stderr)
    sys.exit(1)
print(f"  biopython {Bio.__version__}")
PY

missing=0
for tool in mafft FastTree fasttree augur; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "  ok  $tool -> $(command -v "$tool")"
    else
        # FastTree binary name differs by install (FastTree vs fasttree)
        if [ "$tool" = "FastTree" ] || [ "$tool" = "fasttree" ]; then
            continue
        fi
        echo "MISSING binary: $tool" >&2
        missing=1
    fi
done
if ! command -v FastTree >/dev/null 2>&1 && ! command -v fasttree >/dev/null 2>&1; then
    echo "MISSING binary: FastTree/fasttree" >&2
    missing=1
fi

if [ "$missing" -ne 0 ]; then
    echo >&2
    echo "Install the data-pipeline env:" >&2
    echo "  conda env create -f scripts/panviral/environment.yml" >&2
    echo "  conda activate treesbm-data" >&2
    exit 1
fi

echo "all data-pipeline dependencies found"
