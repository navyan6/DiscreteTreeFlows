"""
Model-independent metric library for the TreeSBM benchmark suite.

Operates on `TreeState` objects and sequence sets; no GPU / model dependency,
so it is unit-testable in isolation.

Submodules:
  trees          tree topology / shape metrics + newick/ete3 I/O
  distributions  distributional distances (Wasserstein, MMD, energy, split-KL)
  sequences      sequence recovery (best-of-K, mutation P/R/F1, coverage@K)
  matching       leaf matching + patristic / sequence-topology coupling
"""
