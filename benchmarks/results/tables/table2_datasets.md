# Table 2 — Datasets and evaluation tasks (no HIV)

**Counted 2026-08-11** from local group CSVs + Betty completed trees (`*rooted.nwk`).  
**# trees / # roots** = trees that finished FastTree + TreeTime + ASR (what train/eval actually load).  
**# sequences** = **leaves** in those completed trees (not internal ASR nodes).  
One reconstructed tree root per group-tree. Eval jobs may subsample *internal* nodes as generation roots (Table 5: 5 from test `group_001`; H3N2 Table 2/4: up to 100).

---

## Paste-ready Table 2

| Dataset | Domain | Root definition | # roots | # trees | # sequences (leaves) | Evaluation split **as implemented** |
|---|---|---|---:|---:|---:|---|
| SARS-CoV-2 Spike | viral evolution | TreeTime ancestor of a single-country, date-contiguous chunk (~300 seqs) | 404 | 404 | 104,260 | **Geographic:** test = **Brazil** (40 trees / 11,907 leaves, 2020-02-28…2025-09-09); val = **Australia** (29 / 8,471). **Not** a future time window. |
| Influenza HA H3N2 | viral evolution | TreeTime ancestor of a date-contiguous collection-year chunk (~400 seqs) | 295 | 295 | 118,000 | **Temporal / future calendar year:** train ≤**2022** (222 / 88,800); val **2023** (28 / 11,200); test **2024** (45 / 18,000). Collection year, not Oct–Sep flu season. |
| Influenza HA H1N1 | viral evolution | TreeTime ancestor of a **location-unit** chunk (≤4-year span, ~50–300 seqs) | 301 | 301 | 55,755 | **Geographic location hold-out** (`data/h1n1`, `prepare_h1n1_geo.py`): 28 held-out GISAID location/country units; test 40 trees / 4,997 leaves (2005–2025). **Not** future-season. A separate temporal 2025 set exists but was **not** used for `h1n1_v2_lit_hotspot` / Table 5. |
| Antibody lineages | affinity maturation | Intended germline/naive; **actual** `root_source=topological` (Zenodo PCP is `no-naive`) | 82 | 82 | 377 | **Held-out Rodriguez clone families** (28 donors). Train for Ab TreeSBM (if used) is separate: Tang+VanWinkle **720 / 80** trees (6,057 / 690 leaves), Rodriguez excluded. Primary 82-family rollout used a **pathogen** ckpt, not Ab-trained. |

HIV Env: skipped.

---

## Train / val / test breakdown (completed trees)

| Dataset | Split | Input groups (fasta/csv) | Completed trees | Leaves in those trees | Years |
|---|---|---:|---:|---:|---|
| COVID geo `data/covid` | train | 396 | **335** usable (336 nwk; group_051 missing `anc_aa`) | 83,882 in 336 nwk | 2019–2026 |
| | val (Australia) | 34 | 29 | 8,471 | 2020–2024 |
| | test (**Brazil**) | 40 | **40** | **11,907** | 2020–2025 |
| H3N2 temporal `data/h3n2` | train ≤2022 | 248 | **222** | 88,800 | 2014–2022 |
| | val 2023 | 30 | **28** | 11,200 | 2023 |
| | test 2024 | 60 | **45** | 18,000 | 2024-01-03…2024-12-02 |
| H1N1 **geo** `data/h1n1` (Table 5 / `h1n1_v2_lit_hotspot`) | train | 218 | 217 | 45,951 | 1998–2026 |
| | val | 45 | 44 | 4,807 | 2004–2026 |
| | test | 40 | **40** | **4,997** | 2005–2025 |
| H1N1 **temporal** `data/h1n1_temporal` (unused for Table 5) | train ≤2023 | 224 | 224 | 37,076 | 2009–2023 |
| | val 2024 | 29 | 27 | 6,328 | 2024 |
| | test 2025 | 33 | 33 | 7,844 | 2025 |
| Antibody eval | Rodriguez held-out | — | **82** | **377** | — |
| Antibody train export | Tang+VanWinkle | — | 720 train / 80 val | 6,057 / 690 | — |

COVID grouped (including failed trees) matches `prepare_covid_geo.py` docstring: val 9,971 / test 11,907 / train 101,330.

---

## H1N1 geo — how location was obtained

**Caption sentence:**  
Location is **not** in per-group `meta.csv` (only `name,date`); it is read from **raw FASTA headers** `ACCESSION \| description \| LOCATION \| DATE \| COUNTRY \| LENGTH`, then whole location units (≥100 seqs) or else country are greedily assigned 80/10/10 so **no place appears in two splits**. Table 5 `group_001` is **Canada: Ontario**, 133 HA, 2009-04-24…2013-12-10.

Not inferred from `A/California/...` strain names. Not GISAID sidecar TSV. Not a temporal split mislabeled as geo.

Raw dumps: `data/h1n1/train/{africa,asia,europe,north_america,oceania,south_america}_h1n1.fasta`.  
Writer strips location: `prepare_h1n1_geo.py` `_emit()` writes `>ACC,DATE` + CSV `name,date` only. **No `SPLIT_PROTOCOL.json`** for geo (unlike temporal).

Test is **28 mixed location units** (Ontario, Iran, pooled Canada, several US states, Finland, Hong Kong, …), **not** a single-country holdout like Brazil.

`h1n1_v2_lit_hotspot` train script hardcodes `data/h1n1/{train,val,test}` (geo). Table 5 coverage job 7539948 used `data/h1n1/test` group_001.

---

## Wording mismatches vs paper screenshot

| Paper phrase | What we actually ran |
|---|---|
| COVID “future time window” | **Brazil geo holdout** (`prepare_covid_geo.py`). Temporal COVID is `data/covid_temporal` (train≤2022 / val 2023 / test 2024–2025) and was **not** Table 5. |
| H3N2 “future season” | **True** as collection-year 2024 vs train≤2022. Not a Northern-Hemisphere Oct–Sep season label. |
| H1N1 “future season” | **False for the used run.** Table 5 / lit-hotspot = **geo location hold-out**. Temporal 2025 (`data/h1n1_temporal`) exists but was not that ckpt/eval. |
| Ab “germline/naive” + “held-out clone family” | Eval **is** 82 Rodriguez families. Roots are **topological** (naive edges already removed upstream). Ab TreeSBM train 720/80 did **not** produce the published 82-family metrics (pathogen ckpt). |

---

## Protocol / audit pointers

| Dataset | Prep / protocol |
|---|---|
| H1N1 geo | `scripts/prepare_h1n1_geo.py` — **no JSON**; `benchmarks/SPLITS.md` §Geographic |
| H1N1 temporal | `data/h1n1_temporal/SPLIT_PROTOCOL.json`, `scripts/prepare_h1n1_temporal.py` |
| H3N2 temporal | `scripts/prepare_h3n2_temporal.py` (default group-size 400, min-year 2014). `SPLIT_PROTOCOL.json` **not** present on Betty/`data/h3n2` |
| COVID geo | `scripts/prepare_covid_geo.py` (`VAL_COUNTRIES={Australia}`, `TEST_COUNTRIES={Brazil}`) |
| COVID temporal | `data/covid_temporal/SPLIT_PROTOCOL.json` (unused for Table 5) |
| Antibody eval | `antibody_benchmark/data/processed/heldout_summary.json`, `FILTER_SUMMARY.md` |
| Antibody train | Betty `data/ab_dasm_trees/export_meta.json` (720/80, `held_out_excluded=rodriguez`) |
| Which run used what | `benchmarks/results/tables/TABLE_METADATA.md`, `benchmarks/SPLITS.md` |
