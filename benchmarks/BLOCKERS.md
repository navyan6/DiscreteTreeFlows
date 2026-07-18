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

## BHV flow (PhylaFlow) — excluded
- Repo: https://github.com/yashaektefaie/PhylaFlow
- Flow matching to a **posterior over trees for a fixed observed alignment** in
  Billera–Holmes–Vogtmann space. No root conditioning, no forward generation;
  it needs the observed sequences and a fixed labeled leaf set.
- **Verdict:** cannot be filled honestly → not in the forward-generation table.

If shown at all, both belong only in a clearly separated *inference-mode
context* note, never in the forward-generation comparison. This matches the
paper's own framing that these methods "cannot perform this task out of the box."

## Non-scientific caveats (still filled, but labeled)
- **TreeSBM** row: its native sampler does not take (N, H); we add an honest
  (N,H) conditioning adapter (`benchmarks/methods/treesbm.py`).
- **ARTreeFormer / PhyloVAE**: topology prior + shared adapters — labeled
  `-adapted`; sequence quality comes from the shared adapter, not the topology
  model.
- **Quartet distance** needs the `tqdist` package; if absent the column is NaN
  (recorded as missing dependency, never faked).
