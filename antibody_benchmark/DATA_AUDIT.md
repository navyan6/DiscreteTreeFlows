# Antibody Affinity-Maturation Benchmark — DATA AUDIT

**Date:** 2026-08-09  
**Repo root:** `/Users/navyanori/Documents/GitHub/DiscreteTreeFlows`  
**Audit paths:** `antibody_benchmark/data/raw/`  
**Betty canonical DASM:** `~/antibody_benchmark_raw/dasm/` (do not redownload if present)

This audit precedes / tracks the fixed-topology recursive rollout benchmark (TreeSBM vs CoSiNE vs DASM+Thrifty vs Thrifty).

---

## 1. Downloads attempted / paths

| Resource | Official location | Local status | Path |
|---|---|---|---|
| CoSiNE HF dataset | `hf://datasets/thematrixmaster/cosine` | **Downloaded** | `antibody_benchmark/data/raw/cosine/hf_dataset/` |
| CoSiNE checkpoint | `hf://thematrixmaster/cosine` → `cosine_dasm.ckpt` (1.7 GB) | **Downloaded** | `antibody_benchmark/data/raw/cosine/checkpoints/cosine_dasm.ckpt` |
| CoSiNE code | `github.com/songlab-cal/cosine` | **Cloned** | `antibody_benchmark/data/raw/repos/cosine/` |
| DASM Zenodo data | Zenodo `17322891` / `dasm-experiments-data.tar.gz` | **Available (Betty canonical + local copy)** | See §1.1 |
| DASM code | `github.com/matsengrp/dasm-experiments` | **Cloned** | `antibody_benchmark/data/raw/repos/dasm-experiments/` |
| netam | `github.com/matsengrp/netam` | **Cloned** | `antibody_benchmark/data/raw/repos/netam/` |
| Thrifty experiments | `github.com/matsengrp/thrifty-experiments-1` | **Cloned** | `antibody_benchmark/data/raw/repos/thrifty-experiments-1/` |
| Thrifty Dryad | DOI `10.5061/dryad.np5hqc044` (`for_dryad.v2.zip`) | **DEFERRED — optional secondary benchmark** | Metadata only: `antibody_benchmark/data/raw/thrifty/dryad_meta.json` |

### 1.1 DASM paths + checksums (Betty canonical)

Configurable via env `AB_DASM_DIR` (or `ANTIBODY_DASM_DIR`) pointing at the extracted `dasm-experiments-data` directory. Config key: `paths.raw_dasm`.

| Artifact | Betty path | size (bytes) | md5 |
|---|---|---:|---|
| Archive | `/vast/home/n/nnori/antibody_benchmark_raw/dasm/dasm-experiments-data.tar.gz` | 128 695 567 | `df38a3488cc8d116bce1250ea18c351c` |
| Archive sha256 | (same file) | — | `18f8710f30af23c479e274277b1db71cdab19087b2e60871f55df9a2f1135265` |
| Extracted root | `/vast/home/n/nnori/antibody_benchmark_raw/dasm/extracted/dasm-experiments-data` | — | — |
| Rodriguez PCP (held-out) | `…/v3/rodriguez-airr-seq-race-prod-NoWinCheck_igh_pcp_2024-11-12_MASKED_NI_noN_no-naive.csv.gz` | 1 770 898 | `05719ad05f8e1e8b3bf99a8daf9b76ff` |

Local Mac mirror (same md5/size; obtained via Betty, not Mac Zenodo):  
`antibody_benchmark/data/raw/dasm/dasm-experiments-data.tar.gz` and `…/extracted/dasm-experiments-data/`.

**Do not redownload** if the Betty archive md5 matches above.

Mac Zenodo direct download may 403 — **not a blocker** (Betty is canonical).

---

## 2. CoSiNE Hugging Face release

### Files

```
train/v1tangCC.txt                 539442 transitions
train/v1vanwinkleheavyTrainCC.txt    99449 transitions
train/v1vanwinklelightTrainCC1m.txt 806639 transitions
train/v1jaffePairedCC.txt           186384 transitions
val/v1tangCC.txt                    112457 transitions
val/v1vanwinkleheavyTrainCC.txt      25536 transitions
val/v1vanwinklelightTrainCC1m.txt   193361 transitions
val/v1jaffePairedCC.txt              23215 transitions
test/v1rodriguezCC.txt               33834 transitions
```

Header line is a count (`N transitions`); each subsequent line is whitespace-separated:

```text
<parent_AA_or_paired> <child_AA_or_paired> <branch_length>
```

Paired Jaffe rows concatenate heavy/light with `.` (e.g. `HEAVY.LIGHT`).

### Columns / recoverable IDs

| Field | Present? |
|---|---|
| parent sequence | yes (AA) |
| child sequence | yes (AA) |
| branch length | yes (float) |
| donor / sample_id | **no** |
| family ID | **no** |
| parent/child node IDs | **no** |
| root/naive flag | **no** |
| leaf flag | **no** |
| V/J gene / CDR coords | **no** |

### Sources (from HF README)

Jaffe 2022, Tang 2022, Vergani 2017, Engelbrecht 2025, Rodriguez 2023.  
Processed with partis + IQ-TREE ancestral reconstruction (same family of PCP construction as DASM).

### Train/val/test

Edge-level split by source file, **not** family-labeled in the text release.

Family counts are **not recoverable** from this release alone (IDs discarded).

### Tree reconstruction

**Not possible** from CoSiNE HF text alone. Edges cannot be linked into clonal trees without `sample_id` / `family` / `parent_name` / `child_name`.

### Alphabet / pairing / BL units

- **AA** (not NT). Light-chain VanWinkle files include leading `X` masks.
- Paired heavy/light retained only in Jaffe files (`.` join).
- Branch length: expected mutations **per site** (CoSiNE README: “calibrated to the expected number of mutations per site”).

---

## 3. DASM Zenodo PCP tables (primary structure source)

### Extracted PCP files (`v3/`)

| File | Source | Rows (incl. header) | Role |
|---|---|---:|---|
| `tang-deepshm-prod-…DXSMVALID.csv.gz` | Tang / DeepSHM productive | 651 936 | train (heavy) |
| `v3convert_vanwinkle-170-igh_…csv.gz` | Engelbrecht/VanWinkle heavy | 124 986 | train (heavy) |
| `wyatt-10x-…_HL.csv.gz` | Jaffe / 10x paired | 209 600 | train (paired HL) |
| `v3convert_vanwinkle-170-igk_…CONCAT…igl….csv.gz` | VanWinkle light | 1 000 001 | train (light) |
| `rodriguez-airr-seq-race-prod-…csv.gz` | Rodriguez RACE-seq | 21 755 | **held-out / perplexity test** |

All published v3 PCP filenames include `no-naive`: **upstream** naive→child edges were already removed in the Zenodo release (not by our pipeline).

### Columns (heavy example)

`sample_id, family, parent_name, parent_heavy, child_name, child_heavy, branch_length, depth, distance, v_gene_heavy, cdr{1,2,3}_codon_{start,end}_heavy, parent_is_naive, child_is_leaf, [j_gene_heavy]`

### Alphabet / pairing / BL

- Sequences are **nucleotide** (ACGTN / masked).
- Heavy/light pairing retained in Jaffe HL file; Tang/Rodriguez/VanWinkle-H are heavy-only.
- `branch_length` and `distance` are continuous floats; consistent with netam / CoSiNE per-site expected-mutation calibration.
- CDR coordinates are **codon indices** on the NT alignment.

### Multi-root cause (Rodriguez held-out)

| Observation | Value |
|---|---:|
| Raw families | 7 769 |
| `parent_is_naive=True` edges in file | **0** |
| Multi-root forests | 578 |
| Of those, identical sequences at all topo-roots | 243 |
| Single-root after resolution | 7 191 |
| Final eligible (≥4 leaves, depth≥2, productive root) | **82** |

**Cause:** Zenodo PCPs are pre-filtered `no-naive`. Removing germline/naive→child edges disconnects sibling lineages that only met at the naive, producing multi-root forests. Identical root sequences across components are consistent with a shared stripped naive.

**Root policy (benchmark):**
1. Prefer an explicit inferred naive/germline (`parent_is_naive`) when present — **do not strip it**.
2. Else, if a unique topological ancestral root remains, retain it for root-conditioned rollout (`root_source=topological`).
3. If multi-root / no usable single root: **exclude** the family. No fabricated edges, no similarity-based ancestry repair.

On Rodriguez Zenodo, all 82 eligible families use `root_source=topological` (naive flags absent in file). Legitimate naive retention is supported in code for any future with-naive PCP; it is not recoverable from this Zenodo release.

### Frozen primary artifacts

| Artifact | Path |
|---|---|
| Trees (JSONL) | `antibody_benchmark/data/processed/benchmark_trees.jsonl` |
| Family IDs | `antibody_benchmark/data/processed/benchmark_family_ids.txt` |
| Attrition CSV | `antibody_benchmark/data/processed/data_audit.csv` |
| Filter summary | `antibody_benchmark/data/processed/FILTER_SUMMARY.md` |

Each JSONL record includes `root_id`, `root_sequence`, `edges`, `branch_lengths`, `leaf_ids`, and eval-only `true_sequences`. Rollout API takes `(model, root_sequence, edges, branch_lengths, seed)` and never receives non-root truths.

### `validate_tree` criteria

1. exactly one root  
2. connected  
3. acyclic  
4. every non-root has indegree 1  
5. every BL > 0  
6. parent/child equal aligned lengths / sequences present  

---

## 4. Thrifty Dryad (neutral SHM secondary)

**Status: DEFERRED — optional secondary benchmark**

- Official archive: `for_dryad.v2.zip` (sha256 `19ddcb2e…`) + `README.md`.
- Intended content: **out-of-frame** and **synonymous-only** PCPs for neutral SHM — **not** the primary affinity-maturation productive benchmark.
- Download blocked here (Dryad API 401). Does **not** block the primary Rodriguez productive rollout benchmark.

---

## 5. Existing TreeSBM Ab work in this repo (do not conflate)

- `data/Homo_sapiens.fasta`, `results/ab_t5/`, scripts such as `scripts/ab_build_clone_trees.py` use OAS / T5-style clone clustering + MAFFT/FastTree.
- That pipeline is **not** the published DASM/CoSiNE PCP tree source. Reuse ANARCI/CDR helpers only; **do not** mix those trees into this benchmark’s held-out set.

---

## 6. Incompatibilities among methods

| Issue | Implication |
|---|---|
| CoSiNE HF has AA + BL but no family/node IDs | Cannot build rollout trees from CoSiNE text; must use DASM CSV structure |
| DASM/Thrifty native alphabet is **NT**; CoSiNE checkpoint is **AA** | Translate NT↔AA at the adapter boundary; metrics A–E primarily on AA Hamming after translation |
| CoSiNE train text is edge-split without family labels | Prefer DASM source files for **family-level** isolation; use Rodriguez as published held-out |
| Upstream `no-naive` multi-root forests | Exclude; do not invent edges |
| Paired HL vs heavy-only | Start primary apples-to-apples on **heavy-only**; keep paired as a separate subset |
| TreeSBM jointly generates topology | For this benchmark **force observed topology+BL**; tree generation is a later experiment |
| CoSiNE Guided Gillespie | **Forbidden** for main benchmark; use unguided `generate_with_gillespie` |
| Thrifty Dryad = neutral / OOF / synonymous | **DEFERRED — optional secondary benchmark** |
| Likelihood comparability | Transition NLL only for models with true normalized p(child\|parent,t); TreeSBM may be N/A |

---

## 7. Recommendation (primary fixed-topology recursive rollout)

**Primary tree source:** DASM Zenodo productive PCP CSVs (`v3/*`), reconstructed into `AntibodyTree` objects with the validation checks above.

**Held-out families:** Rodriguez `rodriguez-airr-seq-race-prod-…csv.gz`, **82** frozen single-root productive heavy trees with ≥4 leaves and depth≥2.

**Alphabet for shared metrics:** translate NT→AA (standard codon table; stop/ambiguous → exclude or mask consistently).

**Secondary (separate):** Thrifty Dryad OOF/synonymous — **DEFERRED**.

**Do not use** CoSiNE HF text as the tree topology source.

---

## 8. Checkpoint / model assets

| Model | Asset | Status |
|---|---|---|
| CoSiNE | `cosine_dasm.ckpt` | local + sync to Betty as needed |
| Thrifty | `netam.pretrained.load("ThriftyHumV0.2-59")` etc. | download-on-demand via netam |
| DASM | `netam.pretrained.load("DASMHumV1.0-4M")` + thrifty neutral | download-on-demand via netam |
| TreeSBM | project checkpoint(s) under `checkpoints/` | configure in YAML |

---

## 9. Blockers / smoke status

1. **Dryad 401** — **DEFERRED — optional secondary benchmark** (neutral SHM only).  
2. **Local Zenodo 403** on this Mac — **non-blocking**; Betty canonical DASM already present (checksums in §1.1).  
3. **Multi-root after upstream no-naive** — understood; 578 families excluded; **82 frozen**. No fabricated edges. Naive retention is implemented; Rodriguez Zenodo has 0 `parent_is_naive` flags → all 82 use `root_source=topological`.  
4. **Betty one-tree smokes (2026-08-10):** Thrifty **PASS**, DASM+Thrifty **PASS**, TreeSBM **PASS** (`~/DiscreteTreeFlows/checkpoints/best.pt`), CoSiNE **PASS** (full ckpt `LABHOME/antibody_benchmark/cosine_ckpts/cosine_dasm.ckpt` md5 `42b14ca6…`; `flash-attn` 2.8.3.post1 in `treesbm` conda). Reports: `antibody_benchmark/results/smoke_reports/`.  
5. Mac CoSiNE/TreeSBM still blocked locally (no CUDA flash-attn / no local `checkpoints/best.pt`).
