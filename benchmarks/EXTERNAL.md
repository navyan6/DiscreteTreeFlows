# External baselines — clone, pin, isolate, sample

The adapted rows (ARTreeFormer, PhyloVAE) use the **official implementations**,
cloned outside DiscreteTreeFlows (on the project volume, not home — see below)
and run in their **own conda envs**. We never re-implement their architectures.
Each is retrained on TRAIN-set size-N topologies and sampled to produce a pool
of newick topologies; the benchmark adapter
(`benchmarks/methods/topology_prior.py`) then consumes that pool + shared
branch-length / sequence adapters. Rows are labeled
`-adapted (topology prior + shared adapters)`.

## Repos + pins

| Method | Repo | Env | Commit |
|---|---|---|---|
| ARTreeFormer | https://github.com/tyuxie/ARTreeFormer | its own `conda env` (`environment.yaml`) | `c2886f49ed8568bfef3d5c058c68ce657e48f3e2` |
| PhyloVAE | https://github.com/tyuxie/PhyloVAE | its own `conda env` (CPU ok, `environment.yaml`) | `7d2867f4e640cc906df39e3c3079ea76d47ab4f0` |

```bash
mkdir -p /vast/projects/pranam/lab/nnori/baselines   # NOT home -- see quota note below
cd /vast/projects/pranam/lab/nnori/baselines
git clone https://github.com/tyuxie/ARTreeFormer
git clone https://github.com/tyuxie/PhyloVAE
conda env create -f ARTreeFormer/environment.yaml
conda env create -f PhyloVAE/environment.yaml
```
Home has a hard 50GB quota (`parcc_quota.py` to check) that's been hit once
already this project from checkpoints/embeddings piling up — clone/train these
somewhere with real headroom.

## The taxon-identity problem (why the pool has to be anonymized)

Both repos hard-require a **single fixed taxon set shared across an entire
training dataset** — they build per-taxon identity embeddings from one `taxa`
list in their own `process_data()` (PhyloVAE: one-hot leaf features via
`torch.eye(ntips)`; ARTreeFormer: `namenum(tree, taxa)`). That's the right
assumption for their original use case (DS1-8: many MCMC tree samples, all
over the *same* fixed species set — posterior uncertainty over one topology).
Our pool is the opposite: many different subtrees from many different H3N2
trees, each with a different set of real leaf strain-ids.

Fix: `export_train_topologies.py` anonymizes every example's leaves to a fixed
generic alphabet `"0".."N-1"` before export. Every example for a given N then
shares the exact same taxon set, so what these models learn is the
exchangeable distribution over N-leaf tree *shapes* — exactly the
unconditional topology prior we want, and consistent with how
`TopologyPriorMethod` only ever consumes bare topology shape anyway (leaf
identity is thrown away and sequence-matched downstream regardless).

## Step 1 — export the training pool (in DiscreteTreeFlows, treesbm env)

```bash
python benchmarks/heldout/export_train_topologies.py \
    --data-dir data/h3n2/train --N 16 32 64 --per-tree 20 \
    --out-dir benchmarks/external_pools
```
Writes `train_topologies_N{16,32,64}.nwk` (bare newick, reference) and
`train_topologies_N{16,32,64}.trprobs` (NEXUS trees block, uniform-weighted —
the exact format `Bio.Phylo.parse(..., 'nexus')` / each repo's
`mcmc_treeprob()` expects). Validated locally by round-tripping a test file
through the identical `Phylo.parse` → `Phylo.write` → `ete3.Tree` path both
repos use. **One thing not in either README:** the weight/rootedness
annotation must come *after* the `=` in the tree line
(`tree t1 = [&W 0.5] [&U] (...);`), not before — `Bio.Nexus` raises a syntax
error otherwise.

## Step 2 — copy the sampling adapters into each repo

```bash
cp benchmarks/external_adapters/artreeformer_sample.py \
   /vast/projects/pranam/lab/nnori/baselines/ARTreeFormer/TDE/
cp benchmarks/external_adapters/phylovae_sample.py \
   /vast/projects/pranam/lab/nnori/baselines/PhyloVAE/
```
These are thin drivers — they load a trained checkpoint and call the repo's
*own* generative method, nothing reimplemented:
- **ARTreeFormer**: `VBPIbase.sample_trees(K)` — true autoregressive sampling,
  builds `K` ete3 trees directly, each attachment conditioned on the partial
  tree so far.
- **PhyloVAE**: *not* autoregressive — `LVMTree.forward()` only evaluates
  existing tree vectors, it never exposes unconditional sampling itself. Real
  recipe, traced from the model code: `z ~ N(0,I)` → `decoder(z)` →
  `cond_prob_mat(...)` gives `[K, ntips-3, 2*ntips-4]` categorical
  distributions, each of the `ntips-3` attachment decisions is drawn
  *independently* given `z` (the continuous latent is what's supposed to
  capture inter-decision correlation), decoded via `vec2tree()` — the exact
  inverse of `tree2vec()`, which their own `process_data()` uses to build
  training vectors. Verified the `edge_mask` valid-position arithmetic matches
  `vec2tree`'s construction step-for-step (row `i` has `2i+3` valid positions,
  matching an `(i+3)`-taxon tree's edge count; total valid positions sums to
  exactly the decoder's flat output dim `(ntips-3)(ntips-1)`) — not yet
  executed against a real trained checkpoint, so sanity-check the first
  sampled pool (right leaf count, no crashes, plausible-looking trees) before
  trusting it at scale.

## Step 3 — train + sample (on the cluster, each repo's own env)

```bash
sbatch scripts/slurm_artreeformer.sh   # GPU, ~12h budget, all 3 N in sequence
sbatch scripts/slurm_phylovae.sh       # CPU (repo's own design), ~12h budget
```
Each: places the `.trprobs` at the path the repo's `process_data()` expects,
runs `process_data`, trains via the repo's own `main.py base.mode=train`, then
samples 300 topologies per N via the copied-in adapter script. Output lands at
`benchmarks/external_pools/sampled/{artreeformer,phylovae}_N{16,32,64}.nwk`.

## Consuming the pool

Already wired into `benchmarks/run_table.py`: if
`benchmarks/external_pools/sampled/*.nwk` exist, `build_methods()` fits a
`BranchLengthAdapter` on `--train-data` and adds both adapted rows
automatically (each N missing a pool is just skipped for that row, not an
error) — no extra flags needed, just run `slurm_artreeformer.sh`/
`slurm_phylovae.sh` before the next `run_table.py` pass. Sequence adapter is
shared JTT (`evolve_pyvolve`, same `--empirical-model` as the native JTT+BD
row) for a like-for-like comparison against the native methods.
