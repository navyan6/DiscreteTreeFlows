# Newick + FASTA

Each tree is a sibling pair with the same stem:

- `*.nwk` — topology (node IDs only; no sequences)
- `*.fasta` — AA sequences; headers must match Newick labels

**Formed (train/val/test):** `data/<dataset>/{train,val,test}/group_*_rooted.nwk` + `group_*_anc_aa.fasta`

**Generated (TreeSBM):** `data/generated/<virus>/{matched,screen,seeds}/group_*_generated*.{nwk,fasta}`

Optional leaf PLL sidecar (same stem): `*.pll.json` — keys are FASTA headers (`node_id|leaf`); Newick tips use `node_id` only.

## Pan-viral data pipeline

Builds formed trees for every eukaryotic virus with enough NCBI genomes
(inventory → CDS extract → MAFFT / FastTree / augur → translate).
Outputs land under `data/panviral/<virus>/{train,test}/`.

### Setup (once)

Do **not** point at someone else's conda path. Create a local env:

```bash
conda env create -f scripts/panviral/environment.yml
conda activate treesbm-data
bash scripts/panviral/check_deps.sh
```

That installs `biopython`, `mafft`, `FastTree`, and `nextstrain-augur`.
Torch/ESM are not required for this pipeline.

Pip-only fallback (still need the three binaries on `PATH`):

```bash
pip install -r scripts/panviral/requirements.txt
```

### Run (Betty / any SLURM login node)

```bash
conda activate treesbm-data
export TREESBM_ROOT=$PWD
export TREESBM_PY=$(which python)
bash scripts/panviral/kickoff.sh
```

One command queues all three stages with SLURM dependencies. Resume is free
(skips finished inventory counts, viruses with a `manifest.json`, and splits
that already have rooted trees).

More detail (knobs, stages, frame checks): `scripts/panviral/README.md`.
