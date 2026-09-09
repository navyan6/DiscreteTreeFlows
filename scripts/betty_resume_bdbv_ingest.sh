#!/bin/bash
# Resume BDBV/filo ingest after NCBI download. Run on Betty login node.
set -euo pipefail
export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:${PATH:-}
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python
cd ~/DiscreteTreeFlows

bash scripts/betty_ingest_filo.sh
