# Pan-viral data pipeline

Goal: not "HIV + SARS-CoV-2 + flu", but every eukaryotic virus with enough
sequence data to build a tree, run through the same align / FastTree / root /
ASR path the COVID datasets already use.

## Kickoff (one command)

On a Betty login node, from the repo root:

```bash
bash scripts/panviral/kickoff.sh
```

That submits stage 1 → 2 → 3 with SLURM dependencies. Later stages wait on
earlier ones, so it is fire-and-forget. Resume is free: inventory caches
counts, stage 2 skips viruses that already have a `manifest.json`, stage 3
skips splits that already have rooted trees.

| Env var | Default | Meaning |
|---|---|---|
| `TREESBM_ROOT` | `$HOME/DiscreteTreeFlows` | repo checkout |
| `TREESBM_PY` | Betty `treesbm` conda python | interpreter with biopython / augur |
| `MIN_COUNT` | `150` | genomes required for a virus to qualify |
| `CHAIN_STAGE3` | `1` | set `0` to stop after the CDS pull |
| `NCBI_API_KEY_FILE` | `~/.ncbi_api_key` | presence only; raises NCBI rate limits |

## Stages

| Stage | What | Scripts |
|---|---|---|
| 1 | Count eukaryotic virus genomes at NCBI | `build_virus_inventory.py`, `slurm_inventory.sh` |
| 2 | RefSeq + antigenic CDS extract, country×time groups | `make_worklist.py`, `fetch_virus_dataset.py`, `slurm_stage2_launcher.sh`, `slurm_fetch_array.sh` |
| 3 | MAFFT / FastTree / augur refine / ancestral / translate | `make_stage3_worklist.py`, `slurm_stage3_launcher.sh`, `slurm_stage3_array.sh`, `scripts/run_all_groups.py` |

Translation mode for new panviral CDS windows is `aligned` with `cds_start=0`
(see `data/TRANSLATION_MODES.md`). Do not reuse flu/HIV `ungapped` here.

## Environment findings (verified 2026-09-08)

These were the open questions that determined the design.

| Question | Answer |
|---|---|
| Does Betty have outbound internet? | Yes — login **and** compute nodes reach NCBI (HTTP 200 from `epyc-3-1`) |
| Tooling present? | `mafft`, `FastTree`, `augur` 33.4.1, `seqkit` all in the `treesbm` conda env |
| Storage headroom? | 19.94 TB of 40 TB used on `/vast/projects/pranam/lab/nnori` |
| SLURM quirks | `--mem` must be a multiple of 5632M; CPU partition is `genoa-std-mem`, not `cpu` |
| NCBI rate limit | 3 req/s unauthenticated (we use 2.5). An API key raises this to 10/s |

## Why not the NCBI Datasets download endpoint

The obvious route — `virus/taxon/{taxid}/genome/download` with
`include_annotation_type=CDS_FASTA` — silently ignores the parameter and returns
a report-only zip. `GENOME_FASTA`, `PROT_FASTA` and `CDS_FASTA` all produce a
byte-identical 305 KB archive containing no sequences. The `dataset_report`
endpoint is also missing the fields we need: across 300 Zika records, zero
carried `isolate.collection_date`, `location` or `host`.

So Datasets is used for one thing only — `total_count` per taxon, which it does
well and cheaply — and everything else goes through Entrez, where GenBank
flatfiles reliably carry `/collection_date`, `/geo_loc_name`, `/host` and the
CDS features.

## Stage 1 — inventory (running, job 8225107)

`scripts/panviral/build_virus_inventory.py`

Downloads the NCBI taxdump once (79 MB) and builds the virus subtree offline,
which avoids ~53k taxonomy API calls. Of 291,212 virus taxa, 53,184 are
species-rank and eukaryote-infecting.

Host range is not in the taxonomy, so eukaryotic status is decided by excluding
clades that are unambiguously prokaryote- or archaea-infecting
(`Caudoviricetes`, `Leviviricetes`, `Microviridae`, `Adnaviria`, the archaeal
families, the phage branches of `Monodnaviria`/`Varidnaviria`) plus any name
containing "phage". The bias is deliberately toward over-inclusion: a stray
phage wastes one tree, a wrongly excluded virus loses a dataset.

Counting is one `total_count` call per species, throttled and appended to a
JSONL cache, so the job is resumable. ~6 h unauthenticated.

Note that `dataset_report` counts include descendants, so species-rank counting
correctly captures Influenza A (11320) and SARS-CoV-2 (under 694009) without
special-casing.

## Stage 2 — reference and protein of interest

The hard part is "ensure proteins of interest are present", because most
GenBank virus records carry **no CDS annotation at all** — of the first three
Zika complete genomes fetched, only one had any features beyond `source`.

Annotating each record is not an option, but we do not have to: every virus
species with data has an annotated RefSeq, confirmed across a deliberately
diverse probe.

| Virus | RefSeq | CDS | Target protein |
|---|---|---|---|
| Zika | `NC_012532.1` | 1 | polyprotein → needs `mat_peptide` |
| Dengue | `NC_001475.2` | 1 | polyprotein → needs `mat_peptide` |
| Measles | `NC_001498.1` | 8 | fusion protein / hemagglutinin |
| Ebola Zaire | `NC_002549.1` | 9 | spike glycoprotein |
| HSV-1 | `NC_001806.2` | 77 | envelope glycoproteins |

So the reference, not the individual record, defines the coordinates — exactly
what `covid_extract_spike.py` already does for Spike. Per virus:

1. Fetch the RefSeq record, read its CDS features.
2. Pick the protein of interest by keyword priority (spike / surface
   glycoprotein / hemagglutinin / envelope / fusion / capsid), falling back to
   the longest CDS when nothing matches.
3. For the two polyprotein families above, descend into `mat_peptide` so we get
   the envelope protein rather than a 10 kb polyprotein.
4. Align every genome to the reference and slice the CDS window by **reference**
   coordinates.

Step 4 is where the Spike frame bug came from, so the slice must stay in
alignment coordinates and translation must go through the fixed gap-aware
`nt_to_aa`.

### Stage 2 validation (2026-09-08)

`scripts/panviral/fetch_virus_dataset.py`, tested on four deliberately awkward
viruses. All resolve the right protein, and QC pass rates run 86–92%.

| Virus | Reference chosen | Target | Length |
|---|---|---|---|
| Ebola Zaire | `NC_002549.1` | spike glycoprotein | 676 aa |
| Zika | `NC_035889.1` | envelope protein E | 504 aa |
| Dengue | `NC_001475.2` | envelope protein E | 493 aa |
| Measles | `NC_001498.1` | hemagglutinin protein | 618 aa |

Zika is the reason `pick_reference` scores several RefSeq candidates instead of
taking the first hit: `NC_012532.1` has a lone polyprotein CDS and no
`mat_peptide` at all, so it yields a useless 3,420 aa target, while
`NC_035889.1` annotates the mature peptides and gives the 504 aa envelope
protein.

Two defects that only showed up under test:

**esearch returns newest-first.** Asking for the first N accessions is a recency
sample, not a sample of the virus. Measles capped at 600 came back with 527
post-2024 records and 36 before — an artefact of ordering that would have
silently destroyed any split by collection year. Fixed by pulling the full UID
list (cheap) and striding through it evenly.

**A fixed 2024 cutoff does not generalise.** It suits continuously sequenced
viruses but empties the test split for outbreak-driven ones: Ebola had 1,036
train against 5 test. The cutoff is now chosen per virus by holding out the most
recent `--holdout-frac` of sequences, which stays strictly forward in time while
adapting to each virus's sampling history.

**The cutoff must land on a window boundary.** Groups are two-year bins but the
cutoff is a single year, so an unsnapped cutoff cuts the bin it lands in: at
cutoff 2018 the 2018–2019 bin sent its 2018 sequences to train and its 2019
sequences to test. For an outbreak still running across that boundary those are
close relatives, which quietly undermines the claim that the test trees are
unseen. `snap_to_window` moves the cutoff down to the end of the previous bin,
which puts the whole contested bin in test — the conservative direction. For
Ebola this moved the cutoff 2018 → 2017 and grew the DRC 2018–2019 test tree
from 80 leaves to 119. Invariants are pinned in
`scripts/panviral/test_cutoff_snap.py`.

The resulting Ebola split is worth noting, because it arrives at the structure
`benchmarks/FILOVIRUS_L_DATA_PLAN.md` argued for on biological grounds — train
on Sierra Leone 2014–2015 (261 leaves, the West Africa outbreak) and forecast
Democratic Republic of the Congo 2018–2019 (80 leaves, Kivu).

## Stage 3 — trees

Unchanged from the existing per-group path: MAFFT the CDS set, FastTree,
`augur refine` for rooting, `augur ancestral` for ASR, then translate with the
gap-aware translator.

## Open scope decisions

- Minimum sequences per virus (currently 150) — sets how many viruses qualify.
- Cap per virus. SARS-CoV-2 alone has millions; trees need thousands, not
  millions, so a per-virus cap sampled across the date range keeps the pull
  bounded.
- An NCBI API key would cut stage 1 from ~6 h to ~1.6 h and speed the much
  larger stage 2 fetch by the same factor.
