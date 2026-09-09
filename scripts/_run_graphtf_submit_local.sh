#!/bin/bash
set -euo pipefail
OUT=/Users/navyanori/Documents/GitHub/DiscreteTreeFlows/scripts/_graphtf_submit_result.txt
{
  echo "START $(date)"
  bash /Users/navyanori/Documents/GitHub/DiscreteTreeFlows/scripts/_paste_on_betty_graphtf.sh
  echo "PASTE_EXIT:$?"
  ssh -o ControlMaster=auto nnori@login.betty.parcc.upenn.edu 'squeue -u nnori -o "%.18i %.12P %.10q %.16j %.8T %R"'
  echo "DONE $(date)"
} >"$OUT" 2>&1
