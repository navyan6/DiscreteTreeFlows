# Table metadata dump (2026-08-11)

Authoritative audit of Tables **2**, **3/4**, **5**, **8**, plus Ab benchmark notes.  
Do not invent numbers; cite artifacts/jobs below.

---

## ε definition (shared)

| Column / name | Code | Meaning |
|---|---|---|
| `coverage` / `eps_frac` | `coverage_at_k(..., eps_frac=0.02)` | **Fractional** Hamming: thresh = `0.02 × L` (ε=2% of length) |
| `coverage_obs_e{E}` | `coverage_at_e(..., e=E)` | **Absolute** AA Hamming ≤ E |

**Table 4 paper Coverage** = `coverage_obs_e5` (ε=**5** absolute Hamming), job **7459412**.  
Not ε=2. The same CSV also has `coverage_obs_e2` for comparison.

**Table 5 morning jobs** (`7539947`/`7539948`) wrote both; paper draft incorrectly used fractional `coverage` (ε=2%).  
Paper-ready Coverage should be **`coverage_obs_e2`** (ε=2 absolute) unless caption says otherwise.

---

## Table 2 — Empirical tree generation baselines (H3N2)

| Field | Value |
|---|---|
| Pathogen | H3N2 HA |
| Split | **Temporal** (`prepare_h3n2_temporal.py`): train ≤2022 / val=2023 / test=**2024** |
| Data | `data/h3n2/{train,val,test}` |
| Years | train 2014–2022 (222 groups); val 2023 (28); test **2024 only** (45 groups, dates 2024-01-03…2024-12-02) |
| N | 16, 32, 64 |
| K / M | K=20 gen, M=20 refs |
| n_roots | max_roots=100 per N (held-out internal roots from test trees) |
| ε | N/A (RF / Quartet / Branch W1 / terminal edit / Tree-KL / Split-KL) |
| TreeSBM ckpt | **`checkpoints/h3n2_v2/best.pt`** (ep57, val=2.3202) |
| Job | **7442467** (CANCELLED time limit 24h; partial CSV pulled) |
| Artifacts | `results/baseline_pull_7442467/results_baselines.csv`, `benchmarks/results/results_baselines.csv`, `table_empirical_*.csv` |
| Script | `scripts/slurm_baselines.sh` → `benchmarks/run_table.py` |
| Coverage roots (first 5 @ N=16, seed=0) | all from **group_002** (2024-01-03…2024-01-05): NODE_0000000, 235, 200, 223, 239 |

Season label: **2024 Northern-Hemisphere calendar-year test window** (collection year 2024), not a flu-season Oct–Sep cut. Train is pre-2023 seasons through 2022.

---

## Table 3 / 4 — H3N2 future-lineage coverage forecasting

| Field | Value |
|---|---|
| Pathogen | H3N2 HA |
| Split / data | Same temporal `data/h3n2/test` (2024) as Table 2 |
| N | **16** |
| K | 10…100 step 10 |
| n_roots | **5** (seed=0; all from **group_002**, dates 2024-01-03…01-05) |
| Legacy ε | `eps_frac=0.02` (fractional; saturates ~1.0) |
| **Paper Coverage (T4)** | **`coverage_obs_e5`** — absolute Hamming ≤ **5** |
| Also computed | `coverage_obs_e{0,1,2,3,5,8,10}` |
| TreeSBM ckpt | **`checkpoints/h3n2_v3_lit_hotspot/best.pt`** (ep57, val=5.1519) |
| Job | **7459412** COMPLETED |
| Script | `scripts/slurm_coverage_curves.sh` → `benchmarks/coverage_curves.py` |
| Artifact | `benchmarks/results/coverage_curves_h3n2_N16_eabs.csv` |
| Log | `benchmarks/results/coverage_curves_7459412.log` |
| Derived tables | `table4_eps5.{md,csv}`, `table4_coverage_vs_eps.csv`, … |

### T4 Coverage@ε=5 (absolute) — job 7459412

| Method | Cov@10 | Cov@100 | mean_min_edit@100 |
|---|---:|---:|---:|
| NeutralBD | 0.800 | 0.800 | 2.65 |
| pLM | 0.800 | 0.800 | 2.65 |
| AR | 0.800 | 0.812 | 2.40 |
| TreeSBM | 0.800 | 0.812 | 2.24 |

Companion **ε=2 absolute** (`coverage_obs_e2` @ K=100): Neutral 0.662 / pLM 0.675 / AR 0.675 / TreeSBM **0.725**.

---

## Table 5 — Viral forecasting (COVID Spike + H1N1 HA)

### Shared protocol

| Field | COVID | H1N1 |
|---|---|---|
| Split | **Geographic** (`prepare_covid_geo.py`): test=**Brazil** | **Geographic** (`prepare_h1n1_geo.py`): held-out locations |
| Data | `data/covid/{train,val,test}` | `data/h1n1/{train,val,test}` |
| Test groups | 40 (`group_001`…`040`) | 40 |
| Test date span | 2020-02-28 … 2025-09-09 | 2005-02-23 … 2025-XX-XX |
| N | 16 | 16 |
| K | 10…100 | 10…100 |
| Coverage script | `scripts/slurm_table5_coverage.sh` | same |
| Enrichment | `scripts/slurm_table5_baselines.sh` / `eval_evescape_enrichment.py` | same |
| Enrichment n_trees | first **20** test trees (groups 1–20) | first **20** |
| L (`max_seq_len`) | **1280** | **566** |

### Morning coverage jobs (fractional paper column — flawed)

| Job | Virus | State | TreeSBM ckpt | Notes |
|---|---|---|---|---|
| **7539947** | COVID | cancelled mid-TreeSBM (~7.7h) | **`covid_v7_pmc_hotspot`** (WRONG vs primary) | Roots = **group_001 ×5** identical leaves → edit=0, cov=1.0 degenerate. pLM OOB: ESM hardcoded L=566 |
| **7539948** | H1N1 | COMPLETED | **`h1n1_v2_lit_hotspot`** ✓ | Roots = **group_001 ×5** (2009-04-24…2013-12-10). OK diversity |

Artifacts: `coverage_curves_{covid,h1n1}_N16_table5.csv`.

### Morning enrichment / ablation companion jobs

| Job | Role | Ckpt |
|---|---|---|
| 7539949 | T5 COVID pLM+AR enrichment | (baselines; no TreeSBM ckpt) |
| 7539950 | T5 H1 pLM+AR enrichment | — |
| 7539951–54, 7539980, 7546797 | Table 8 gen ablations | **`covid_v5_mutrec`** |
| 7539981 | Table 8 retrain λ_br=0 | trains `covid_ablate_no_bl` |
| 7546809 | COVID Tree-KL / Split-KL / W1 / pLM NLL | **`covid_v5_mutrec`** — **FAILED** RateHeads 1280 vs default 566 |

### Best checkpoints (enrichment @ mrs=0.5)

| Virus | Best ckpt | Evidence |
|---|---|---|
| COVID | **`checkpoints/covid_v5_mutrec/best.pt`** | mut≈0.069–0.080, aa\|hit≈0.34–0.38; **v7** mut≈0.045, aa\|hit≈0.157 (worse) |
| H1N1 | **`checkpoints/h1n1_v2_lit_hotspot/best.pt`** | used correctly in 7539948 |

**Resubmit needed for COVID coverage:** morning used **v7** + degenerate group_001 + ESM L=566 bug.

### Resubmitted e-abs jobs (this audit)

| Job | Name | Purpose |
|---|---|---|
| **7575015** | `t5_cov_covid_eabs` | COVID Coverage e-abs; ckpt=**v5_mutrec**; ESM+TreeSBM L=1280; diverse one-root-per-group filter; out=`coverage_curves_covid_N16_table5_eabs.csv` |
| **7575016** | `t5_cov_h1_eabs` | H1 re-run with same e-abs protocol / out=`…_h1n1_…_eabs.csv` (ckpt still `h1n1_v2_lit_hotspot`) |
| **7575017** | `t8_base_covid` | COVID baselines Tree-KL/Split-KL/W1/pLM NLL with `--max-seq-len 1280` |

Code fixes synced to Betty: `run_table.py --max-seq-len`, `coverage_curves.py` ESM(`max_len=…`), diverse-root filter, `slurm_baselines.sh` / `slurm_table5_coverage.sh`.

### Old vs new Coverage numbers

#### H1N1 — job 7539948 (existing; absolute already in CSV)

| Method | Cov@10 frac | Cov@100 frac | **Cov@10 e=2** | **Cov@100 e=2** | min_edit@100 |
|---|---:|---:|---:|---:|---:|
| pLM | 0.975 | 0.975 | **0.362** | **0.362** | 4.04 |
| AR | 0.975 | 0.975 | **0.362** | **0.462** | 3.88 |
| TreeSBM | 0.988 | 1.000 | **0.375** | **0.525** | 3.43 |

Use **e=2** columns for paper; fractional ε=2% is nearly saturated.

#### COVID — job 7539947 (DO NOT USE)

| Method | Cov@10/100 frac & e2 | min_edit | Issue |
|---|---|---|---|
| pLM / AR / Neutral | 1.0 | 0.0 | identical GT leaves in group_001 |
| TreeSBM | not flushed | — | job cancelled |

Await **7575015** for paper-ready COVID Cov@10/100 @ e=2.

### Table 5 coverage roots (morning jobs)

**COVID 7539947** (all group_001, Brazil, dates ~2020-02-28…2020-XX-XX, L=1273, mean_edit=0):

- NODE_0000074, NODE_0000177, NODE_0000218, NODE_0000223, NODE_0000188

**H1N1 7539948** (all group_001, dates 2009-04-24…2013-12-10, L=566):

- NODE_0000076, NODE_0000095, NODE_0000117, NODE_0000067, NODE_0000019

### Table 5 enrichment groups (max_trees=20)

COVID Brazil test groups **1–20** (date spans 2020-02-28 → 2021-09-14; see Betty meta CSVs).  
H1N1 geo test groups **1–20** (mixed locations/years 2005–2025; group sizes 52–244).

---

## Table 8 — COVID ablations

| Field | Value |
|---|---|
| Pathogen | SARS-CoV-2 Spike |
| Ckpt (all eval-only rows) | **`checkpoints/covid_v5_mutrec/best.pt`** |
| Data | `data/covid/test` (+ train for baselines) |
| Mask / EVEscape | `results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt`, `data/covid/evescape_spike_rbd.pt` |
| mrs / n_steps / max_trees | 0.5 / 100 / 20 |
| Script | `scripts/slurm_table8_ablations.sh` |
| Job map | `benchmarks/results/tables/table8_job_ids.json` |

| Job | Ablation |
|---|---|
| 7539951 | `--ablate-bridge` |
| 7539952 | `--ablate-tree-context` |
| 7539953 | poisson_ref branching |
| 7539954 | `--ablate-branch-length-head` |
| 7539980 | full |
| 7546797 | `--ablate-internal-node-seqs` |
| 7539981 | retrain `LAMBDA_BR=0` → `covid_ablate_no_bl` |
| 7546809 → **7575017** | Tree-KL / Split-KL / Branch W1 / pLM NLL |

Coverage@100 per ablation: not run (heavy). Enrichment fills mut recall / aa\|hit / cons.

---

## Antibody benchmark (brief)

| Field | Value |
|---|---|
| Data | 82 families, N=20 (`configs/full.yaml`) |
| Jobs | thrifty 7517208, dasm 7517209, treesbm 7517210, cosine 7517211, eval 7517212 |
| TreeSBM caveat | 7517210 used pathogen `checkpoints/best.pt`, **not** Ab-trained |
| Ab retrain | 7539619 → `checkpoints/ab_dasm_v1/best.pt` (loss collapsed); rollout 7539827 |

---

## Metric definitions (quick)

| Metric | Definition |
|---|---|
| Coverage@K (frac) | Fraction of GT leaves with min Hamming ≤ `eps_frac·L` to some gen leaf in first K trees |
| Coverage@K (e-abs) | Same with absolute threshold e |
| mean_min_edit | Mean over GT leaves of min Hamming to gen pool |
| site_recall / mut_recovery / aa\|hit / cons | From `mutation_pr_f1` / `positional_recovery` |
| Tree-KL / Split-KL | Topology distribution KL (`benchmarks/metrics/distributions.py`) |
| Branch W1 | Wasserstein on branch lengths (`branch_w_all`) |
| pLM NLL | From baselines / enrichment path when wired |

---

## Code edits applied (Betty synced)

1. `benchmarks/run_table.py`: `--max-seq-len` → ESM + `TreeSBMMethod`
2. `benchmarks/coverage_curves.py`: `ESM(..., max_len=args.max_seq_len)`; diverse one-root-per-group selection
3. `scripts/slurm_baselines.sh`: `MAX_SEQ_LEN` env → CLI
4. `scripts/slurm_table5_coverage.sh`: default COVID ckpt → **v5_mutrec**; eabs out paths

---

## Refresh note (2026-08-15 ~06:05 ET Betty pull)

- **Table 7 / repo T7:** ESM-C **7627569/70** COMPLETED (n=20); ProGen2 **7627918** COMPLETED (n=20; supersedes **7627571** n=3). Additive cons/aa|hit/antigenic filled. Cov@100 for ESM-C/ProGen2 **7628260/61/62 RUNNING** — see `table7_refproc.md` / PASTE §2.
- **HIV C.1:** baselines **7626609/10**, cov **7626614/15** COMPLETED. Tables under `tables/hiv_{geo,temporal}/`; TreeSBM topo n_roots geo=8 / temporal=6.
- **Ab OAS:** train **7627289** FAILED after writing `best.pt`; Rod.82 **7627290** DependencyNeverSatisfied. No new Track C samples.
