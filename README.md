# TreeSBM: Tree-Valued Schrödinger Bridge Matching

A generative model for influenza phylogenetic trees. Given a root ancestral sequence, TreeSBM simulates viral evolution by growing a bifurcating phylogeny whose leaves are descendant sequences with biologically plausible mutations and branching topology.

## Motivation

Phylogenetic trees are the central object in molecular evolution — they encode how viruses spread, mutate, and compete over time. Existing tools (FastTree, RAxML) *infer* trees from observed sequences but cannot *generate* new ones. TreeSBM learns a generative model directly from empirical flu trees, enabling simulation of counterfactual evolutionary trajectories from any ancestral sequence.

The model is trained on 253 time-calibrated influenza trees spanning 7 subtypes (H3N2, H1N1 HA, H1N1 NA, swine H3N2, avian H1N1, flu B Yamagata, flu B Victoria) built from GISAID sequences.

## Method

TreeSBM uses **Schrödinger Bridge Matching** to learn a stochastic process that starts at a root-only tree (T=0) and ends at a full observed phylogeny (T=1). At each intermediate time t, a partial tree T_t is sampled by time-cutting the observed tree, and the model is supervised on the rates that drive T_t toward T1.

The learned rate function decomposes as:

```
log R_θ = log R0 + c_θ
```

where R0 is ESM-2's masked language model head (a biologically grounded prior over amino acid substitutions) and c_θ is a learned correction conditioned on tree context and bridge time t.

**Architecture:**
- **NodeEncoder**: fuses per-node ESM-2 embedding (320-dim, mean-pooled) + structural features (depth, children, is_leaf) + Laplacian positional encoding → 128-dim
- **TreeEncoder**: 4-layer graph transformer with temporal causal attention mask (nodes attend only to earlier ancestors) and branch-length-weighted edges, conditioned on bridge time t
- **RateHeads**: per active leaf, predicts mutation logits [L×20] (conditioned on per-position ESM-2 log-probs), branching rate λ, branch length, and stop probability

**Training losses:**
| Loss | Signal |
|------|--------|
| L_seq | Cross-entropy on log R_θ_mut vs. T1 amino acids, weighted by 1/(1-t) (Doob h-transform) |
| L_top | Poisson NLL on branching rate vs. log(total T1 descendants + 1) |
| L_br | MSE on branch length vs. mean T1 child branch lengths |
| L_stop | BCE on stop probability vs. whether node is a T1 terminal leaf |
| L_pll | ESM-2 pseudo-log-likelihood regularizer on current T_t sequences |

## Data Pipeline

Raw FASTA (GISAID) → date-stratified split into groups of ~400 sequences → MAFFT alignment → FastTree → augur/TreeTime (root + time calibration) → ancestral sequence reconstruction → ESM-2 embeddings + reference rates precomputed per node.

## Setup

```bash
conda env create -f environment.yml
conda activate treesbm
```

Requires MAFFT, FastTree, and augur (available via the nextstrain conda channel).

## Usage

**Precompute embeddings** (run once per dataset):
```bash
python scripts/precompute_plm.py --data data/train
python scripts/precompute_ref_rates.py --data data/train
```

**Train:**
```bash
sbatch scripts/slurm_train.sh
# or locally:
python scripts/train.py --data data/train --epochs 300 --patience 50
```

**Generate a tree from a root sequence:**
```bash
python scripts/generate_tree.py \
    --checkpoint checkpoints/best.pt \
    --root-seq <amino_acid_sequence> \
    --n-steps 50 \
    --branch-rate-scale 6.0 \
    --max-seq-len 566 \
    --output generated_tree.nwk
```

**Evaluate on held-out test trees:**
```bash
python scripts/eval_test_set.py \
    --checkpoint checkpoints/best.pt \
    --data data/train \
    --split checkpoints/split_indices.json \
    --n-steps 50 \
    --max-leaves 200 \
    --branch-rate-scale 6.0
```

## Repository Structure

```
src/
  bridge/          # SampleBridgeState (Algorithm 2) and training losses
  treeencoder/     # NodeEncoder, ESM-2 embedder, Laplacian PE, structural features
  networks.py      # TreeEncoder and RateHeads
  dataset.py       # TreeDataset
  tree_state.py    # TreeState: mutable tree representation
scripts/
  train.py         # Training loop (Algorithm 1)
  generate_tree.py # Controlled generation (Algorithm 4)
  eval_test_set.py # Batch evaluation on held-out test set
  precompute_plm.py
  precompute_ref_rates.py
  slurm_train.sh
  slurm_eval.sh
data/train/        # Per-group FASTA, Newick, metadata, PLM/ref-rate caches
checkpoints/       # best.pt + split_indices.json
```
