# External baselines — clone, pin, isolate, sample

The adapted rows (ARTreeFormer, PhyloVAE, PhylaFlow) use the **official
implementations**, cloned outside DiscreteTreeFlows (on the project volume, not
home — see below) and run in their **own conda envs**. We never re-implement
their architectures. Each is retrained on **our H3N2 TRAIN-set** size-N
topologies (or, for PhylaFlow, on H3N2 training alignments / posterior contexts
that yield size-N topology samples — **not** PhylaFlow’s paper DS1–8) and
sampled to produce a pool of newick topologies; the benchmark adapter
(`benchmarks/methods/topology_prior.py` / `phylaflow_native.py`) then consumes
that pool + shared branch-length / sequence adapters. Adapted rows are labeled
`-adapted (topology prior + shared adapters)`.

**PhylaFlow is NOT native forward generation** — it is posterior-basin transport
in BHV space for a fixed alignment. The `phylaflow_adapted` row uses only the
sampled topology pool + shared adapters, same honesty pattern as ARTreeFormer /
PhyloVAE. Native `phylaflow` keeps PhylaFlow BLs (rescaled to H) + shared JTT.

## Repos + pins

| Method | Repo | Env | Commit |
|---|---|---|---|
| ARTreeFormer | https://github.com/tyuxie/ARTreeFormer | its own `conda env` (`environment.yaml`) | `c2886f49ed8568bfef3d5c058c68ce657e48f3e2` |
| PhyloVAE | https://github.com/tyuxie/PhyloVAE | its own `conda env` (CPU ok, `environment.yaml`) | `7d2867f4e640cc906df39e3c3079ea76d47ab4f0` |
| PhylaFlow | https://github.com/yashaektefaie/PhylaFlow | its own `conda env` (see repo; GPU typical) | pin at clone time |

```bash
export LABHOME=/vast/projects/pranam/lab/nnori
mkdir -p $LABHOME/baselines   # NOT home -- see quota note below
cd $LABHOME/baselines
git clone https://github.com/tyuxie/ARTreeFormer
git clone https://github.com/tyuxie/PhyloVAE
git clone https://github.com/yashaektefaie/PhylaFlow
conda env create -f ARTreeFormer/environment.yaml
conda env create -f PhyloVAE/environment.yaml
# PhylaFlow: follow its README / environment spec at clone time
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

Submit **one N per job** (16 / 32 / 64). Packing all three into one job risks
starving N=64 of walltime (observed on PhyloVAE). Jobs for different N can run
concurrently.

```bash
# GPU — ARTreeFormer (own `artreeformer` conda env)
sbatch scripts/slurm_artreeformer.sh 16
sbatch scripts/slurm_artreeformer.sh 32
sbatch scripts/slurm_artreeformer.sh 64

# CPU — PhyloVAE (own `phylovae` conda env; repo is CPU-only by design)
sbatch scripts/slurm_phylovae.sh 16
sbatch scripts/slurm_phylovae.sh 32
sbatch scripts/slurm_phylovae.sh 64

# GPU — PhylaFlow (own env; H3N2 bank train/sample first — NOT DS1–8)
# Default QOS in script header is mig-max; override with sbatch --qos=...
# Convert dumps → pools only after H3N2-trained dumps exist (see PhylaFlow section).
sbatch --qos=mig-max scripts/slurm_phylaflow.sh 16
sbatch --qos=mig-max scripts/slurm_phylaflow.sh 32
sbatch --qos=mig-max scripts/slurm_phylaflow.sh 64
```
Each ARTreeFormer / PhyloVAE job: places the `.trprobs` at the path the repo's
`process_data()` expects (directory name reuses the DS1–8 *layout* but the
trees are our anonymized H3N2 train topologies), runs `process_data`, trains
via the repo's own `main.py base.mode=train`, then samples 300 topologies via
the copied-in adapter script. Output lands at
`benchmarks/external_pools/sampled/{artreeformer,phylovae}_N{N}.nwk`.

## PhylaFlow — native row (`phylaflow`) + adapted (`phylaflow_adapted`)

PhylaFlow learns hybrid flow matching in BHV tree space to transport random
starting trees toward **posterior basins for a fixed observed alignment**. It
does not accept `(root_seq, N, H)` for forward generation.

### What “DS1–8” is (and why Table 2 does **not** use it)

**DS1–DS8** are the classic fixed-taxon MrBayes posterior banks from PhylaFlow’s
own paper (`short_run_data_DS1-8/`, `golden_run_data_DS1-8/`,
`fixed_path_artifacts/DS{1..8}/`). Their `./launch_ds_local.sh ds1`…`ds8` /
`joint` recipes reproduce *PhylaFlow Table 1 / 12* — not TreeSBM Table 2.

For a fair RF / Quartet / Branch-W1 / Terminal-edit comparison against TreeSBM
and the BD/PLM rows on **H3N2 held-out roots** (`data/h3n2/test`, n≈97),
PhylaFlow must be trained (or finetuned) and sampled on **our H3N2 train
alignments / posterior contexts** that yield size-N trees — the same empirical
pool distribution as ARTreeFormer/PhyloVAE (`export_train_topologies.py` on
`data/h3n2/train`). **Do not** train Table-2 pools with `launch_ds_local.sh ds2`
(or any other DS1–8 case). DS1–8 roots are optional only if you want to
reproduce PhylaFlow’s paper numbers separately.

Same honesty as other adapted methods: still not root-conditioned forward gen;
we consume anonymized sampler trees as a size-N prior + shared JTT sequences.

**Betty clone path (lab volume, not home):**

```bash
export LABHOME=/vast/projects/pranam/lab/nnori
export DTF=~/DiscreteTreeFlows
# One paste after Duo SSH:
bash $DTF/scripts/_paste_betty_phylaflow_native.sh
# When H3N2 dumps/pools ready:
bash $DTF/scripts/_paste_betty_phylaflow_native.sh submit
```

| Row | What it uses from PhylaFlow | Sequences | Branch lengths |
|---|---|---|---|
| **`phylaflow`** (native / Table 2 target) | Official sampler trees (topo + BLs) from **H3N2** bank | shared JTT | PhylaFlow BLs, rescaled to H |
| `phylaflow_adapted` | Same pool topologies | shared JTT | shared `BranchLengthAdapter` |

For the benchmark we still need an offline pool (never invent trees):

1. Build a **custom H3N2 size-N bank** in PhylaFlow’s layout (alignments +
   posterior / path artifacts for train subtrees), under e.g.
   `$PHYLAFLOW_DATA_ROOT/h3n2_N{N}/` — *not* `short_run_data_DS1-8/`.
   Start from anonymized train topologies:
   `python benchmarks/heldout/export_train_topologies.py --data-dir data/h3n2/train …`
   Then follow PhylaFlow’s bank-construction scripts (fixed-path anchors,
   phyla embeddings, metric start table) pointed at those H3N2 cases.
2. Train PhylaFlow on that H3N2 bank (own env, GPU, qos=mig-max) with a
   size≈N config — **not** `./launch_ds_local.sh ds*`.
3. Dump sampled trees with PhylaFlow's
   `scripts/evaluate_per_dataset_sample_kl.py --dump-trees`.
4. Convert dumps → anonymized pool via
   [`benchmarks/external_adapters/phylaflow_sample.py`](external_adapters/phylaflow_sample.py)
   (leaves → `"0".."N-1"`, **keep BLs** by default, one tree per line).
5. Land at `benchmarks/external_pools/sampled/phylaflow_N{16,32,64}.nwk`.

**Do not fake pools.** If no dump/checkpoint exists, leave the row absent
(`run_table.py` skips missing pools).

### Cluster protocol (mig-max)

```bash
export LABHOME=/vast/projects/pranam/lab/nnori
export DTF=/vast/home/n/nnori/DiscreteTreeFlows   # or your checkout

# 0) Clone under lab volume (not home)
git clone https://github.com/yashaektefaie/PhylaFlow $LABHOME/baselines/PhylaFlow
cd $LABHOME/baselines/PhylaFlow
pip install -r requirements.txt
# H3N2 custom bank roots (NOT DS1–8):
export PHYLAFLOW_DATA_ROOT=$LABHOME/baselines/phylaflow_h3n2_data
export PHYLAFLOW_ARTIFACT_ROOT=$LABHOME/baselines/phylaflow_h3n2_artifacts
export PHYLAFLOW_OUTPUT_ROOT=$LABHOME/baselines/PhylaFlow/outputs_h3n2

cp $DTF/benchmarks/external_adapters/phylaflow_sample.py \
   $LABHOME/baselines/PhylaFlow/phylaflow_sample.py

# 1) Anonymized H3N2 train topologies (treesbm env) — same pool as ARTree/PhyloVAE
cd $DTF
python benchmarks/heldout/export_train_topologies.py \
    --data-dir data/h3n2/train --N 16 32 64 --per-tree 20 \
    --out-dir benchmarks/external_pools

# 2) Build PhylaFlow H3N2 bank from those train cases (PhylaFlow scripts;
#    alignments / MrBayes or equivalent posteriors / fixed-path anchors /
#    phyla embeddings). Layout under $PHYLAFLOW_DATA_ROOT/h3n2_N{N}/ …
#    Do NOT copy or train on short_run_data_DS1-8 for Table 2.

# 3) Train PhylaFlow on the H3N2 bank (own env, GPU, qos=mig-max) — one N
cd $LABHOME/baselines/PhylaFlow
# Use / adapt a config with leaf count ≈ N and data paths → h3n2_N{N}
#   sbatch --qos=mig-max --partition=b200-mig45 --gres=gpu:1 --wrap \
#     'cd $LABHOME/baselines/PhylaFlow && python -m run.run --config configs/h3n2_N16.yaml'
# (FORBIDDEN for Table 2: ./launch_ds_local.sh ds1…ds8)

# 4) Sample + dump trees from the H3N2-trained ckpt
python scripts/evaluate_per_dataset_sample_kl.py \
    --config configs/h3n2_N16.yaml \
    --checkpoint $PHYLAFLOW_OUTPUT_ROOT/<ckpt>.ckpt \
    --sample-config <sample_metrics.yaml> \
    --output-dir samples/treesbm_N16 \
    --num-samples 50 --dump-trees

# 5) Convert dumps → benchmark pool (or submit the wrapper job)
python phylaflow_sample.py \
    --input samples/treesbm_N16/tree_dumps \
    --ntips 16 --n-samples 300 \
    --out $DTF/benchmarks/external_pools/sampled/phylaflow_N16.nwk

# Symlink / convert helper (exits 1 with protocol if nothing ready):
cd $DTF
sbatch --qos=mig-max scripts/slurm_phylaflow.sh 16
sbatch --qos=mig-max scripts/slurm_phylaflow.sh 32
sbatch --qos=mig-max scripts/slurm_phylaflow.sh 64
```

Primary Table-2 row: **`phylaflow`** via
`benchmarks/methods/phylaflow_native.py` (`NativePhylaFlowMethod`). Optional
ablation row: **`phylaflow_adapted`** via `TopologyPriorMethod` (shared BL
adapter). Both auto-register in `run_table.py` when
`benchmarks/external_pools/sampled/phylaflow_N{N}.nwk` exists.

## EVEscape (Marks) — escape enrichment evaluator (not a generator)

Official EVEscape scores single-AA mutants as:

```
evescape = log σ(z(fitness_eve)/T_f)
         + log σ(z(accessibility_wcn)/T_a)
         + log σ(z(dissimilarity)/T_d)
```

So the **final** score is a sum of three log-probabilities: typically
**negative and outside [0,1]** (flu H1 / Spike RBD means ≈ −2.3). Do **not**
expect [0,1]; that applies only to intermediate logistic terms. Raw EVE
(`fitness_eve`) is a different column (flu ≈ −11; Spike ≈ −6).

### Preferred path: Marks precomputed CSVs (no PDB needed)

```bash
# From https://github.com/OATML-Markslab/EVEscape/tree/main/results/summaries_with_scores
mkdir -p data/evescape_official
curl -L -o data/evescape_official/flu_h1_evescape.csv \
  https://raw.githubusercontent.com/OATML-Markslab/EVEscape/main/results/summaries_with_scores/flu_h1_evescape.csv
curl -L -o data/evescape_official/spike_rbd_evescape.csv \
  https://raw.githubusercontent.com/OATML-Markslab/EVEscape/main/results/summaries_with_scores/spike_rbd_evescape.csv

# H1N1 HA (L=566) — keep official scale (--no-standardize is default)
python scripts/prepare_evescape.py \
  --csv-path data/evescape_official/flu_h1_evescape.csv \
  --output data/evescape_h1n1_ha.pt \
  --max-seq-len 566 --pathogen flu_h1 \
  --ref-seq-file data/evescape_h1n1_ref_seq.txt   # or --data data/h1n1/train

# COVID Spike RBD on full-spike frame (L=1280)
python scripts/prepare_evescape.py \
  --csv-path data/evescape_official/spike_rbd_evescape.csv \
  --output data/covid/evescape_spike_rbd.pt \
  --max-seq-len 1280 --pathogen covid_spike_rbd \
  --ref-seq-file data/covid/evescape_ref_seq.txt
```

**Never** pass `--standardize` for official `evescape` columns (legacy COVID
tensors did; that zero-centers the paper scale and makes `model_evescape≈0`).

Pathogen separation: flu_h1 CSV → H1 trees only; Spike RBD → COVID only.
There is **no** official H3N2 EVEscape summary in the Marks release.

### Optional: recompute from PDB (only if you need a new structure/WT)

EVEscape repo (`process_protein_data.py` → `evescape_scores.py`) needs:

| Pathogen | PDB (paper) | Chains | EVE indices | WT fasta |
|---|---|---|---|---|
| Flu H1 HA | **1RVX** (`1rvx_no_HETATM.pdb`) | A,B (+ trimer A–F) | I4EPC4 evol indices | `A0A2Z5U3Z0_9INFA.fasta` |
| Spike RBD | 6VXX / 6VYB / 7BNN / 7CAB | see EVEscape | P0DTC2 pre-2020 | `SPIKE_SARS2.fasta` |

You do **not** need to send a PDB for H1 if Marks `flu_h1_evescape.csv` is enough.

### Eval metrics (`scripts/eval_evescape_enrichment.py`)

EVEscape is defined **per mutation** (official CSV log-scale values). We never
multiply mutations into a strain-level product/sum.

| Key | What it averages |
|---|---|
| `evescape_mean_antigenic_muts` / `gt_evescape_mean_antigenic_muts` | Mean EVEscape of root→leaf AA changes **at lit/antigenic mask sites** (`--lit-hotspot-mask`: H1 Sa/Sb/Ca1/Ca2/Cb∪guidance; COVID PMC; H3 lit when used). Preferred flu antigenic KPI. |
| `random_baseline_evescape_antigenic` | Mean of nonzero tensor entries restricted to the same antigenic columns. |
| `evescape_mean_all_scored_muts` / `model_evescape` | Mean over **all** root→leaf muts with nonzero tensor entry (legacy companion; scored region ≈ RBD or full HA). |
| `random_baseline_evescape` | Mean of all nonzero tensor entries (no site mask). |

Also emits `n_*_evescape_antigenic_muts_scored` and mean antigenic-mut count per
leaf. `model_frac_rbd` is a legacy name; `model_frac_scored` is the clearer alias.
Pathogen separation: never pass COVID RBD EVEscape / PMC mask with flu HA `L`.

## EVE (Marks) — sequence / mutation baseline (not a tree generator)

EVE evolutionary indices score single-AA mutants (log-likelihood ratio vs WT).
We use them as an **external evaluator** parallel to EVEscape (Table 7 / B.3),
not as a forward tree method.

### Produce `--eve-scores` `.pt` (offline)

```bash
export LABHOME=/vast/projects/pranam/lab/nnori
git clone https://github.com/OATML-Markslab/EVE $LABHOME/baselines/EVE
cd $LABHOME/baselines/EVE
conda env create -f protein_env.yml && conda activate protein_env

# MSA: download from https://evemodel.org/ (preferred) or build via Jackhmmer/EVcouplings
# Then (see examples/):
bash examples/train_VAE.sh              # → VAE checkpoint for your protein
bash examples/compute_evol_indices.sh   # → CSV with mutations + evol_indices
```

Convert CSV → TreeSBM column frame (`[L, 20]`, AA order `ACDEFGHIKLMNPQRSTVWY`):

| Protein | Typical out | `--data` | `--max-seq-len` | Notes |
|---|---|---|---|---|
| Spike (COVID) | `data/covid/eve_spike.pt` | `data/covid/train` | `1280` | Full spike frame used by TreeSBM; use `--position-offset` if EVE MSA is mature/RBD-only |
| HA (H3N2) | `data/h3n2/eve_ha.pt` | `data/h3n2/train` | `566` | Align EVE WT to group root via `prepare_eve_scores.py` |
| HA (H1N1) | `data/h1n1/eve_ha.pt` | `data/h1n1/train` | `566` | Same converter; check match rate ≥ 0.90 |

```bash
# HA example
python scripts/prepare_eve_scores.py \
    --csv-path $LABHOME/baselines/EVE/results/evol_indices/<HA>_20000_samples.csv \
    --output data/h3n2/eve_ha.pt \
    --data data/h3n2/train --ref-from-group 1 --max-seq-len 566

# Spike example
python scripts/prepare_eve_scores.py \
    --csv-path $LABHOME/baselines/EVE/results/evol_indices/<Spike>_20000_samples.csv \
    --output data/covid/eve_spike.pt \
    --data data/covid/train --ref-from-group 1 --max-seq-len 1280
```

Alternative: download precomputed per-protein EVE tables from
[evemodel.org](http://evemodel.org/) and convert the mutant score table with the
same script (columns `mutations` / `evol_indices` or pass `--score-col`).

### Score TreeSBM generations

```bash
# Local / interactive
python scripts/eval_eve_baseline.py \
    --checkpoint checkpoints/covid_v3_cons/best.pt --data data/covid/test \
    --max-seq-len 1280 --eve-scores data/covid/eve_spike.pt \
    --mutation-rate-scale 0.3 --n-steps 100 --max-trees 20

# Cluster (GPU for generation; EVE lookup is CPU once .pt exists)
sbatch --qos=mig-max scripts/slurm_eval_eve.sh \
    checkpoints/covid_v3_cons/best.pt data/covid/test data/covid/eve_spike.pt 1280
sbatch --qos=mig-max scripts/slurm_eval_eve.sh \
    checkpoints/h3n2_v2/best.pt data/h3n2/test data/h3n2/eve_ha.pt 566
```

Reported JSON keys: `recovered_gt_eve`, `random_aa_eve`,
`recovered_minus_random_eve`, `rate_eve_correlation.{pearson_r,spearman_r}`.

## Consuming the topology pools

Already wired into `benchmarks/run_table.py`: if
`benchmarks/external_pools/sampled/*.nwk` exist, `build_methods()` fits a
`BranchLengthAdapter` on `--train-data` and adds all three adapted rows
automatically (each N missing a pool is just skipped for that row, not an
error) — no extra flags needed, just run `slurm_artreeformer.sh` /
`slurm_phylovae.sh` / `slurm_phylaflow.sh` before the next `run_table.py` pass.
Sequence adapter is shared JTT (`evolve_pyvolve`, same `--empirical-model` as
the native JTT+BD row) for a like-for-like comparison against the native methods.
