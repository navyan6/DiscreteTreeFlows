#!/bin/bash
# ONE paste on Betty after Duo SSH (cwd ~/DiscreteTreeFlows).
# Sync from laptop first if needed, then:
#   bash scripts/_paste_betty_ops_all.sh
set -euo pipefail
cd ~/DiscreteTreeFlows
bash scripts/_paste_betty_siteaa_reeval.sh
bash scripts/_paste_betty_track_a.sh
bash scripts/_paste_betty_graphtf_j1.sh
