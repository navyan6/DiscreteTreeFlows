# Run from laptop (scp + remote submit). Or paste the remote half on login03.
set -euo pipefail
REPO=/Users/navyanori/Documents/GitHub/DiscreteTreeFlows
HOST=nnori@login.betty.parcc.upenn.edu
scp -o ControlMaster=auto \
  "$REPO/scripts/validate_transformer.py" \
  "$REPO/scripts/slurm_validate_transformer.sh" \
  "$REPO/scripts/check_column_entropy.py" \
  "$REPO/scripts/slurm_check_column_entropy.sh" \
  "$REPO/scripts/_cluster_submit_graphtf_probe.sh" \
  "$HOST:~/DiscreteTreeFlows/scripts/"
ssh -o ControlMaster=auto "$HOST" 'bash ~/DiscreteTreeFlows/scripts/_cluster_submit_graphtf_probe.sh'
