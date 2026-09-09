# Rows that cannot be filled honestly in the forward-generation table

The task is **forward generation**: given only (root sequence, N, H), produce a
rooted descendant tree with *novel* terminal sequences, never seeing the real
descendants. Two named methods are phylogenetic **inference** methods that
require the observed terminal sequences and therefore cannot perform this task
without becoming a different method. They are **excluded** from the strict table.

## PhyloGFN — excluded
- Repo: https://github.com/zmy1116/phylogfn
- A GFlowNet whose **reward is the phylogenetic likelihood / parsimony of a tree
  given the alignment's terminal sequences**. Sampling is per-alignment and
  driven entirely by those sequences.
- To run it forward we would have to (a) give it the hidden test sequences
  (leakage — forbidden), or (b) replace its reward with something sequence-free,
  which re-purposes the GFlowNet into a different model. Its per-alignment
  training also does not yield a transferable, sequence-free, size-N topology
  prior.
- **Verdict:** cannot be filled honestly → not in the forward-generation table.

If shown at all, PhyloGFN belongs only in a clearly separated *inference-mode
context* note, never in the forward-generation comparison. This matches the
paper's own framing that these methods "cannot perform this task out of the box."

## Adapted topology prior only (NOT native forward generation)

These methods are phylogenetic **inference / posterior transport** models. They
do **not** perform root-conditioned forward generation out of the box. We use
them honestly as **unconditional size-N topology priors** (retrain or sample on
anonymized train-set topologies, then consume the pool via
`TopologyPriorMethod` + shared branch-length / sequence adapters). Rows are
labeled `-adapted`; sequence quality comes from the shared adapter, not the
topology model.

- **ARTreeFormer** (`artreeformer_adapted`): see `benchmarks/EXTERNAL.md`.
- **PhyloVAE** (`phylovae_adapted`): see `benchmarks/EXTERNAL.md`.
- **PhylaFlow** (`phylaflow` native + optional `phylaflow_adapted`): repo
  https://github.com/yashaektefaie/PhylaFlow — hybrid flow matching to a
  **posterior basin for a fixed observed alignment** in Billera–Holmes–Vogtmann
  space. No root conditioning; it needs observed sequences / bank context.
  **Native Table-2 row** = official sampler topo+BL (rescaled to H) + shared JTT
  sequences — still not root-conditioned forward gen. Train/sample on **our
  H3N2 train bank** (same empirical distribution as Table 2 / ARTreeFormer /
  PhyloVAE), **not** PhylaFlow’s paper DS1–8. Requires clone under
  `$LABHOME/baselines/PhylaFlow`, custom H3N2 bank + dumps →
  `benchmarks/external_pools/sampled/phylaflow_N{N}.nwk` (see
  `scripts/_paste_betty_phylaflow_native.sh`, `EXTERNAL.md`).

## Non-scientific caveats (still filled, but labeled)
- **TreeSBM** row: its native sampler does not take (N, H); we add an honest
  (N,H) conditioning adapter (`benchmarks/methods/treesbm.py`).
- **Quartet distance** needs the `tqdist` package; if absent the column is NaN
  (recorded as missing dependency, never faked).
