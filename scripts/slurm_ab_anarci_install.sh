#!/bin/bash
#SBATCH --job-name=ab_anarci_install
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=8
# omit --mem: PARCC filter requires exact 5632MB/CPU on genoa-std-mem
#SBATCH --time=2:00:00
#SBATCH --output=logs/ab_anarci_install_%j.log
#SBATCH --error=logs/ab_anarci_install_%j.log
#
# Install ANARCI (pip) + HMMER (build from eddylab tarball) into treesbm.
# Betty login often has no conda CLI; pip anarci alone is NOT enough (needs hmmscan).

set -euo pipefail

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
ENV_PREFIX=/vast/home/n/nnori/.conda/envs/treesbm
BIN="$ENV_PREFIX/bin"

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Start: $(date) host=$(hostname)"

if ! $PYTHON -c "import anarci" 2>/dev/null; then
  echo "pip install anarci"
  $PYTHON -m pip install anarci
fi

if ! command -v hmmscan >/dev/null 2>&1; then
  echo "Building HMMER 3.4 into $BIN"
  WORK=${TMPDIR:-/tmp}/hmmer_build_$$
  mkdir -p "$WORK"
  cd "$WORK"
  if [[ ! -f hmmer-3.4.tar.gz ]]; then
    curl -fsSL -o hmmer-3.4.tar.gz http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
  fi
  tar xf hmmer-3.4.tar.gz
  cd hmmer-3.4
  ./configure --prefix="$ENV_PREFIX"
  make -j"${SLURM_CPUS_PER_TASK:-4}"
  make install
  cd ~/DiscreteTreeFlows
  rm -rf "$WORK"
fi

export PATH="$BIN:$PATH"
command -v hmmscan
$PYTHON -c "import anarci; print('ANARCI_OK', anarci.__file__)"

$PYTHON - <<'PY'
from anarci import anarci
seq = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINSSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARPPRNYYDRSGYYQRAEYFQHWGQGTLVTVSS"
r = anarci([("t", seq)], scheme="imgt", output=False, assign_germline=True)
assert r[0][0], "ANARCI returned empty numbering — HMMER/DB issue"
print("ANARCI_RUNTIME_OK")
PY

$PYTHON -c "import dendropy; print('DENDROPY_OK', dendropy.__version__)"
command -v mafft
command -v FastTree
echo "Done: $(date)"
