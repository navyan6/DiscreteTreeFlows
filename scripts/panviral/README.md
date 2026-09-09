# Pan-viral data pipeline

Build trees for every eukaryotic virus with enough NCBI genomes: inventory → fetch & extract the antigenic CDS → align / FastTree / root / ASR.

## Run it

On a Betty login node, from the repo root:

```bash
git pull
bash scripts/panviral/kickoff.sh
```

That one command queues all three stages with SLURM dependencies. You can walk away.

```bash
squeue -u $USER                 # see the jobs
tail -f logs/panviral/*.log     # watch progress
```

Resume is automatic: re-running `kickoff.sh` skips finished inventory counts, viruses that already have a `manifest.json`, and splits that already have rooted trees.

## Optional knobs

| Variable | Default | Meaning |
|---|---|---|
| `MIN_COUNT` | `150` | min genomes for a virus to qualify |
| `CHAIN_STAGE3` | `1` | set `0` to stop after CDS fetch (no trees yet) |
| `TREESBM_ROOT` | `$HOME/DiscreteTreeFlows` | repo path on the cluster |
| `TREESBM_PY` | Betty `treesbm` conda python | python with biopython / mafft / FastTree / augur |

Example — fetch only, require ≥100 genomes:

```bash
MIN_COUNT=100 CHAIN_STAGE3=0 bash scripts/panviral/kickoff.sh
```

Optional NCBI API key (never put it in the repo): write it once to `~/.ncbi_api_key` with `chmod 600`. If that file exists and is non-empty, the launcher raises concurrency; the key is never printed or passed on the command line.

## What each stage does

1. **Inventory** — count NCBI genomes for eukaryotic virus species → `data/panviral/virus_inventory.json`
2. **Fetch** — RefSeq + antigenic protein CDS, grouped by country × time window → `data/panviral/<virus>/{train,test}/`
3. **Trees** — MAFFT, FastTree, `augur refine`, `augur ancestral`, translate → `group_*_rooted.nwk` + `group_*_anc_aa.fasta`

Outputs land under `data/panviral/`. Logs under `logs/panviral/`.

## Before you trust a new virus's proteins

New viruses must use the **aligned** translator with `cds_start=0` (already set in stage 3). If you add a virus outside this pipeline, check the frame first:

```bash
python scripts/panviral/audit_frame_offset.py data/panviral/<virus>/train
python scripts/panviral/probe_cds_start.py data/panviral/<virus>/train
```

Do not copy flu/HIV `ungapped` mode onto reference-anchored CDS windows.
