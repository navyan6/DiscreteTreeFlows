# VOC threat recovery panel (Track A)

Goal: show TreeSBM recovers defining mutations of pandemic VOCs from regional
ancestral roots, with honest bundle / top-k metrics and EVEscape annotations.

**Read first:** [`SPIKE_FRAME_BUG.md`](SPIKE_FRAME_BUG.md) (what was wrong),
[`NUMBERS_THAT_CHANGE.md`](NUMBERS_THAT_CHANGE.md) (what moved),
[`PAPER_NUMBERS_TO_UPDATE.md`](PAPER_NUMBERS_TO_UPDATE.md) (what it means for the
ICLR submission), [`LEAKAGE.md`](LEAKAGE.md) (the standard we apply).

All panel numbers below are **post-fix**. Mutation calling now maps reference
positions through a pairwise alignment (`voc_threat_lib.seq_at`), so
deletion-carrying lineages — Alpha, Beta, Delta, Omicron BA.1 — are read at the
right residue. Anything quoted from before 2026-09-07 is stale.

## The standard

A case supports a recovery claim when the **root sequence was never seen in
training** and the **tree is held out**. Shared mutations between train and test
are expected and fine — a pandemic corpus is full of D614G — and the model is
supposed to learn escape trajectories that reuse common mutations. What it must
not have is the specific ancestor it is asked to forecast from.

Two columns below encode this: `root_seen` (root Spike string appears in train)
and `novel` (bundle-carrying leaves whose Spike is *not* byte-identical to any
training Spike).

## Canonical VOC table

[`benchmarks/voc_threat_panel.json`](../../benchmarks/voc_threat_panel.json) —
signatures, defining bundles, origin countries, date windows, candidate groups.

## Panel A — held-out Brazil (geographic test split)

The strongest cases, because the split was designed for this and the entire
country is absent from training.

| Case | VOC | `bundle_frac` | Notes |
|---|---|---|---|
| **g003** | Gamma | **0.717** | Figure 4/5 case study. Root 1273 aa, 98.7% frame-consistent. P.1 has no Spike deletion, so this case is the control for the whole audit — it did not move. |
| g011 | Delta | **0.418** | Was 0.033 before the fix (12.5×). Caution: only 59.5% of its leaves share the root's frame, so treat column-wise metrics on this tree carefully. |
| g007 | Delta | **0.053** | Was 0.010 (5.3×). 95.3% frame-consistent. |

Manifest: `cases_manifest_clean.json`.

## Panel B — fresh NCBI roots, unseen populations

Pulled by `download_sars2_voc_roots_ncbi.py` with accession-level exclusion
against `seen_accessions.txt`, then built into trees on Betty. Verified by
`check_voc_roots_overlap.py` → `voc_roots_seq_overlap.json`.

| Population | VOC | tree | leaves | bundle leaves | novel | root unseen |
|---|---|---|---|---|---|---|
| **usa_iota** | Iota | g002 | 299 | 79 | **79 (100%)** | **yes** |
| **colombia_mu** | Mu | g001 | 291 | 125 | **125 (100%)** | **yes** |
| **usa_epsilon** | Epsilon | g003 | 300 | 120 | 33 | **yes** |
| **usa_epsilon** | Epsilon | g002 | 300 | 118 | 30 | **yes** |
| **uk_alpha** | Alpha | g003 | 260 | 123 | 20 | **yes** |
| usa_iota | Iota | g003 | 260 | 80 | 80 (100%) | no |
| usa_iota | Iota | g001 | 300 | 69 | 69 (100%) | no |
| south_africa_beta | Beta | g001 | 300 | 143 | 26 | no |
| uk_alpha | Alpha | g002 | 300 | 124 | 26 | no |
| uk_alpha | Alpha | g001 | 298 | 123 | 15 | no |
| usa_epsilon | Epsilon | g001 | 300 | 118 | 27 | no |
| ~~india_delta~~ | Delta | g001 | 150 | 1 | 0 | no |

### What changed with the fix

Alpha and Beta looked unusable and are in fact among the best-populated trees:
Beta g001 went from 3 bundle-carrying leaves to **143**, Alpha g001 from 2 to
**123**. The pre-fix screen was reading their signature positions two residues
off because of the 69-70del / 144del and 241-243del deletions.

Epsilon and Iota are unchanged, as expected — neither lineage carries a Spike
deletion, so they serve as the second control alongside Gamma.

### Submission tiers

**Tier 1 — headline prospective claims.** Unseen root, every bundle-carrying leaf
novel: `usa_iota g002`, `colombia_mu g001`. USA is the cleanest population in the
panel because North America is entirely absent from the training corpus.

**Tier 2 — unseen root, partly novel leaves:** `usa_epsilon g002`, `g003`,
`uk_alpha g003`. Report the novel fraction in the caption rather than burying it.

**Tier 3 — seen root, case-study framing only:** `south_africa_beta g001`,
`uk_alpha g001/g002`, `usa_epsilon g001`, `usa_iota g001/g003`. These carry the
most bundle leaves but the root was in training, so they illustrate rather than
demonstrate.

**Dropped:** `india_delta g001` — one bundle leaf and a seen root.

## Scripts

| Script | Role |
|--------|------|
| `scripts/voc_threat_lib.py` | Mut parsing / **alignment-aware** scoring / EVEscape lookup |
| `scripts/screen_voc_threat_trees.py` | Screen observed trees for VOC content + prospective roots |
| `scripts/eval_voc_threat_recovery.py` | Score generated leaves (exact / acquired / bundle / top-k) |
| `scripts/prepare_voc_origin_cases.py` | Copy case folders + write manifest |
| `scripts/download_sars2_voc_roots_ncbi.py` | Pull unseen VOC-origin genomes from NCBI |
| `scripts/prepare_covid_voc_roots.py` | Group pulled seqs into temporally stratified trees |
| `scripts/check_voc_roots_overlap.py` | Sequence-level leakage check on new roots |
| `scripts/audit_train_leakage.py` | Re-audit VOC carriers in train (alignment-aware) |
| `scripts/check_spike_frame_shift.py` | Direct test of the coordinate-shift hypothesis |
| `scripts/check_tree_frame_consistency.py` | Per-tree leaf-vs-root frame consistency |
| `scripts/audit_frame_all_datasets.py` | Frame exposure across covid / h3n2 / h1n1 / HIV |
| `scripts/slurm_covid_voc_roots_pipeline.sh` | Betty: spike extract → group → align → tree → ASR |
| `scripts/slurm_voc_threat_panel.sh` | Betty: generate via `eval_single_tree` + VOC eval |
| `scripts/betty_voc_threat_go.sh` | One-shot sync + submit |

## Betty submit

```bash
# Panel A (Brazil, held out)
MANIFEST=results/voc_threat_panel/cases_manifest_clean.json \
  CKPT_NAME=covid_v7_pmc_hotspot MRS=0.5 MAX_LEAVES=250 N_STEPS=160 \
  bash scripts/betty_submit_voc_threat_panel.sh

# dry-run first
DRY=1 MANIFEST=results/voc_threat_panel/cases_manifest_clean.json \
  bash scripts/betty_submit_voc_threat_panel.sh
```

Jobs write `cases/<VOC>_gXXX_<split>/generated.{nwk,fasta}` and
`voc_eval_<VOC>.json`.

## Next steps

1. Screen Panel B trees for prospective roots (`screen_voc_threat_trees.py`),
   then build a Tier-1/Tier-2 manifest.
2. Submit Panel A and the Tier-1 Panel B cases for generation + VOC eval.
3. Open: reference-coordinate masks (RBD / PMC literature / EVEscape tensor) are
   applied to padded columns, so they land on the wrong residues for the 77.5% of
   COVID roots that are not 1273 aa. This is the last unquantified piece and it
   touches training, not just reporting — see §5 of `PAPER_NUMBERS_TO_UPDATE.md`.
