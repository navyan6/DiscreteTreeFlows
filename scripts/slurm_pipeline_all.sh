#!/bin/bash
#SBATCH --job-name=treesbm_pipeline
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm_pipeline_%j.log
#SBATCH --error=logs/slurm_pipeline_%j.log

export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
PYTHON=/vast/home/n/nnori/.conda/envs/treesbm/bin/python

cd ~/DiscreteTreeFlows
mkdir -p logs

echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURMD_NODENAME"
echo "Start:  $(date)"
echo "Python: $($PYTHON --version)"

# Run all 7 prefixes in parallel (4 workers each = 28 concurrent groups)
$PYTHON scripts/run_all_groups.py --prefix h3n2_swine_all          --group-offset 48  --workers 4 &
$PYTHON scripts/run_all_groups.py --prefix avian_h1n1_2010_2020_HA --group-offset 55  --workers 2 &
$PYTHON scripts/run_all_groups.py --prefix h1n1_human_ha_2010_2017 --group-offset 56  --workers 4 &
$PYTHON scripts/run_all_groups.py --prefix human_h1n1_NA_2005_2015 --group-offset 106 --workers 4 &
$PYTHON scripts/run_all_groups.py --prefix human_h1n1_2015_2018    --group-offset 156 --workers 4 &
$PYTHON scripts/run_all_groups.py --prefix fluB_yamagata_alltime   --group-offset 206 --workers 4 &
$PYTHON scripts/run_all_groups.py --prefix fluB_victoria_all       --group-offset 238 --workers 4 &

wait
echo "All pipelines done: $(date)"
