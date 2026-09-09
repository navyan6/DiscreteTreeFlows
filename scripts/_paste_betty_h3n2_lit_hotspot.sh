#!/bin/bash
# Betty paste: H3N2 lit-hotspot train (nmicrobiol201658 mask).
# PREPARED ONLY — do not auto-submit. Sync masks + code first.
#
# STRICT SEPARATION: uses results/flu_mutfreq_vs_lit/* only.
# Never point MUT_HOTSPOT_MASK at results/covid_mutfreq_vs_lit/.
#
# Prereqs on Betty:
#   1. rsync results/flu_mutfreq_vs_lit/ (masks + COORD_VALIDATION.json)
#   2. data/h3n2/{train,val,test} with group_*_anc_aa.fasta (pipeline done)
#   3. Optional: recompute masks on Betty:
#        python scripts/validate_ha_coords.py --data data/h3n2/train --max-seq-len 566 \
#            --out-dir results/flu_mutfreq_vs_lit
#        python scripts/compute_flu_mut_freq_vs_lit.py --data data/h3n2/train \
#            --max-seq-len 566 --out-dir results/flu_mutfreq_vs_lit
#
# Then paste / run one of the sbatch lines below (edit qos/partition as needed).

set -euo pipefail
cd ~/DiscreteTreeFlows

FLU_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt
test -f "$FLU_MASK" || { echo "MISSING $FLU_MASK — sync or regenerate first"; exit 1; }

# --- Preferred: curated nmicrobiol lit mask + force (antigenic sites) ---
# Mirrors covid_v7_pmc_hotspot but H3N2 paths / L=566 / flu mask only.
cat <<EOF
# === PASTE ON BETTY (lit hotspot; not submitted by agent) ===
HOTSPOT=1 \\
CKPT_DIR=checkpoints/h3n2_v3_lit_hotspot \\
MUT_HOTSPOT_MASK=results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt \\
MUT_HOTSPOT_FRAC= \\
MUT_HOTSPOT_TOPK= \\
MUT_HOTSPOT_WEIGHT=5 \\
MUT_HOTSPOT_FORCE=1 \\
LAMBDA_MUT=12 LAMBDA_CONS=0.5 LAMBDA_SEMI=0.05 \\
MUT_NORMALIZE=count PSSM_GATE=1 MUT_AA_EMB=1 \\
sbatch --qos=mig-max --job-name=h3n2_lithot \\
  --export=ALL,HOTSPOT,CKPT_DIR,MUT_HOTSPOT_MASK,MUT_HOTSPOT_FRAC,MUT_HOTSPOT_TOPK,MUT_HOTSPOT_WEIGHT,MUT_HOTSPOT_FORCE,LAMBDA_MUT,LAMBDA_CONS,LAMBDA_SEMI,MUT_NORMALIZE,PSSM_GATE,MUT_AA_EMB \\
  scripts/slurm_h3n2_train_lit_hotspot.sh
# === END PASTE ===
EOF

echo ""
echo "If slurm_h3n2_train_lit_hotspot.sh is missing, use train.py directly:"
cat <<'EOF2'
sbatch --qos=mig-max --job-name=h3n2_lithot --gres=gpu:1 --cpus-per-task=8 --mem=64G --time=24:00:00 \
  --wrap='export PATH=/vast/home/n/nnori/.conda/envs/treesbm/bin:$PATH
cd ~/DiscreteTreeFlows
python -u scripts/train.py \
  --data data/h3n2/train --val-data data/h3n2/val --test-data data/h3n2/test \
  --max-seq-len 566 --epochs 100 --patience 50 --bridge-c 1.0 \
  --lambda-mut 12 --lambda-cons 0.5 --mut-normalize count --lambda-semi 0.05 \
  --per-site-pos-emb --use-site-entropy \
  --use-entropy-loss-weighting --use-entropy-cons-weighting \
  --entropy-source empirical --entropy-weight-alpha 3.0 --entropy-weight-alpha-cons 1.0 \
  --pssm-gate --mut-aa-emb \
  --mut-hotspot-mask results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt \
  --mut-hotspot-weight 5 --mut-hotspot-force \
  --ckpt-dir checkpoints/h3n2_v3_lit_hotspot'
EOF2

echo ""
echo "Eval after train (explicit flu mask; never PMC):"
cat <<'EOF3'
python scripts/eval_evescape_enrichment.py \
  --checkpoint checkpoints/h3n2_v3_lit_hotspot/best.pt \
  --data data/h3n2/test --max-seq-len 566 \
  --lit-hotspot-mask results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt \
  --mutation-rate-scale 0.3 \
  --out checkpoints/eval_enrichment_h3n2_v3_lit_hotspot_mrs0.3.json
EOF3
