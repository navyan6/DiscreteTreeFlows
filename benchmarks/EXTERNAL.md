# External baselines — clone, pin, isolate, sample

The adapted rows (ARTreeFormer, PhyloVAE) use the **official implementations**,
cloned into sibling folders and run in their **own conda envs**. We never
re-implement their architectures. Each is retrained on TRAIN-set size-N
topologies and sampled to produce a pool of newick topologies; the benchmark
adapter (`benchmarks/methods/topology_prior.py`) then consumes that pool +
shared branch-length / sequence adapters. Rows are labeled `-adapted (topology
prior + shared adapters)`.

## Repos + pins (record the exact commit you use)

| Method | Repo | Env |
|---|---|---|
| ARTreeFormer | https://github.com/tyuxie/ARTreeFormer | its own `conda env` per its README |
| PhyloVAE | https://github.com/tyuxie/PhyloVAE | its own `conda env` (CPU ok) |

```bash
cd ..                      # sibling of DiscreteTreeFlows
git clone https://github.com/tyuxie/ARTreeFormer && cd ARTreeFormer && git rev-parse HEAD  # PIN THIS
git clone https://github.com/tyuxie/PhyloVAE     && cd PhyloVAE     && git rev-parse HEAD  # PIN THIS
```
Write both commits at the bottom of this file.

## Producing the topology pool (offline, in the repo's own env)

1. **Training topologies (from TreeSBM train trees):** for each N ∈ {16,32,64},
   sample the induced size-N subtree *topologies* from `data/h3n2/train` (reuse
   `benchmarks.heldout.build_examples.induced_subtree`, drop branch lengths +
   sequences → unlabeled newick). Save as the `.trprobs`/newick-set format the
   repo expects (both ARTreeFormer TDE and PhyloVAE consume tree-sample sets).
2. **Train** the method on those topologies, per N, following the repo's TDE
   instructions (its architecture + objective unchanged).
3. **Sample** ≥ K topologies per N from the trained model → newick file
   `pool_artreeformer_N{N}.nwk` (one newick per line). Same for PhyloVAE.

## Consuming the pool in the benchmark
```python
pool = {N: Path(f"pool_artreeformer_N{N}.nwk").read_text().splitlines() for N in (16,32,64)}
method = TopologyPriorMethod("artreeformer_adapted", pool, bl_adapter, seq_adapter_fn)
```
`bl_adapter` = `BranchLengthAdapter().fit(train_trees)`; `seq_adapter_fn` = the
shared sequence adapter (same for both adapted methods, e.g. shared JTT via
`evolve_pyvolve` or shared PLM via `evolve_plm`).

## Pinned commits (fill in)
- ARTreeFormer: `______`
- PhyloVAE: `______`
