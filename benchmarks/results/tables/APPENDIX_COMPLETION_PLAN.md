# Appendix completion plan (ICLR TreeSBM PDF)

**Created:** 2026-08-15  
**PDF:** `TreeSBM__ICLR_ (8).pdf` (29 pp; appendix starts p.19)  
**Deadline window:** **5 days** remaining — do **not** try to fill everything at once.  
**Rule:** never invent numbers; cite artifacts/jobs. Empty cells stay **—** with an explicit blocker.  
**Additive layers (on top of existing metrics, never replace):** `cons_retention`, antigenic-site proportion (`lit_hotspot_mut_frac` / flu/HIV V-loop frac), `aa_acc_given_hit` (aa|hit). Apply on viral enrichment / forecasting / ablation / reference tables where sequence gens exist.

Cross-check sources: `PASTE_TABLES_2026-08-14.md`, `TABLE_METADATA.md`, `table7_refproc.md`, `table8_ablations.md`, `table6_ab_track_c.md`, `antibody_tracks.md`, `table_c3_*`, `hiv_*` tables, `PAPER_TABLES.md`, `EVAL_METRICS_LEAF_VS_TREE.md`, main-text Tables 2–8.

## Standing SOP

On every result pull (cluster job COMPLETED → scp/rsync → table fill), update this plan in the same pass:

- Mark checklist items done, paste key metric numbers, clear or refresh blockers, record job IDs
- Do not leave this plan stale relative to `PASTE_TABLES` / `table_*.md`

---

## 0. Executive inventory (appendix only)

| Appendix table / fig | Paper columns (all) | Status | Priority (5-day) |
|---|---|---|---|
| **A.1** Dataset construction | Root / Leaf / Internal / BL / Filtering | **DONE** (text) | — |
| **A.2** Split definitions | Train / Test / Generalization | **DONE** (text) | — |
| **B.1–B.3** Metric glossaries | Meaning / Direction | **DONE** (text) | — |
| **C.1** Full tree gen by dataset | Tree-KL, Split-KL, RF, Quartet, Branch W1, pLM NLL | **H1N1 swap** HIV Env → Influenza HA H1N1; H1 cells **--** pending **7646375–77**; COVID/H3/Ab intact | 2026-08-17 |
| **C.2** Future-lineage full | Cov@100, Cov@500, Mut recall, Clade recall, Min dist | **NEARLY FULL** — HIV pLM + Ab Cov/Clade/Min filled; COVID TSBM Cov@500 **7635911** RUNNING | Day 2 + 2026-08-16 ✓ |
| **C.3** Horizon | Cov@100, Mut recall, Clade recall, Min dist, Tree-KL | TreeSBM LaTeX [`table_c3_horizon.tex`](table_c3_horizon.tex): clade COVID/H3 **0.000**, HIV geo long **0.097**; empty buckets COVID short / H3 short / HIV geo med; Tree-KL — | Day 1 ✓ |
| **D.1** pLM reference full | Params, Fitness, Tree-KL, Split-KL, Cov@100, Mut, pLM NLL | **ENRICH + Cov DONE** (ESM-C/ProGen2 Cov **7628260–62** filled) | Day 1 ✓ |
| **D.2** Fitness/branching priors | Cov@100, Tree-KL, Branch W1 | **FILLED** (ESM2 W1 **7628576/77** → 5.20e-5 / 5.32e-5) | Day 4 ✓ |
| **E.1** Bridge matching | Tree-KL, Cov@100, Mut recall | **PARTIAL** (term-only / w/o term filled; no_doob **7628590** RUNNING) | Day 2 ✓ / Day 5+ |
| **E.2** Architecture ablations | Tree-KL, Cov@100 | **QUEUED** mut/stop-head trains (see `table_e2_job_ids.json`) | Day 2 ✓ |
| **E.3** Sampling sensitivity | Cov, Diversity, Tree-KL, Runtime | **TEMP/MAXL FILLED** (**7639276–79**); K=500 **7635911** RUNNING | Day 5 ✓ / 2026-08-16 |
| **F.1** Data scaling | Tree-KL, Cov@1000, Mut, pLM NLL | **MISSING** | **Defer** (blocker: no train fraction jobs) |
| **F.2** Cross-domain transfer | Tree-KL, Cov@100, Mut, Branch W1, pLM NLL | **MISSING** | **Defer** |
| **G.1** Likelihood ranking | True / Random / Matched / Shuffled ranks | **MISSING** | **Defer** |
| **Fig G.1** Coverage calibration | figure | **MISSING** | **Defer** |
| **H.1 / H.2** Qualitative galleries | figures | **DOCUMENTED** (tree_viz selection) | Day 5 ✓ |
| **I.1** Compute cost | Params, train time, mem, time/root | **FILLED** (`table_i1_compute.tex`; GPU mem —) | 2026-08-17 ✓ |
| **J.1** Linear-probe embeddings | Leaf/Int … HamPar (9 cols × 6 emb) | **DONE** (H3N2+COVID probes) | Day 1 ✓ |
| **J.2** Hyperparameters | Value column | **DONE** (`table_j2_j3_hyperparams.tex` + `.md`) | Day 5 ✓ / LaTeX **2026-08-17** |
| **J.3** Sampling settings | Value column | **DONE** (`table_j2_j3_hyperparams.tex` + `.md`) | Day 5 ✓ / LaTeX **2026-08-17** |
| **K** Algorithms | pseudo-code | **DONE** (text) | — |

**Main-text tables that appendix must stay consistent with** (not appendix IDs, but paste blockers): Table 5 Ab (paper T5) = Rod.82 Track C; paper Table 6 = repo T7 refproc; paper Table 7 = repo T8 ablations; paper Table 8 panel design = mostly empty.

---

## 1. Day-by-day schedule (5 days)

### Day 1 — Paste what exists + close D.1 Cov gap *(today)* — **DONE except Cov pull**

1. [x] **C.3 paste** from `table_c3_horizon.md` (jobs **7627191–95** COMPLETED) → also in `PASTE_TABLES` §3c. Leave Tree-KL as —.
2. [x] **Ab OAS Rod.82** paste `treesbm_ab_oas` into paper Table 5 / Track C (`table6_ab_track_c.md` + `PASTE_TABLES` §3b; job **7627919** COMPLETED). **No further OAS PLM precompute needed for Rod.82.**
3. [x] **D.1 / paper Table 6:** ProGen2 enrich **7627918** COMPLETED (n=20). Cov@100 submitted (overnight):
   - **7628260** `esmc_nofit` coverage — **RUNNING** (pull Day 2)
   - **7628261** `esmc_fit` coverage — **RUNNING** (pull Day 2)
   - **7628262** `progen2_fit` coverage — **RUNNING** (pull Day 2)
4. [x] Refresh `PASTE_TABLES` / `table7_refproc.md` with ProGen2 n=20 numbers + additive cons / aa|hit / antigenic. Cov cells stay **—** until 7628260–62 finish.
5. [x] **J.1** confirmed FINAL in `PAPER_TABLES.md` / wave4 artifacts (`PASTE_TABLES` §3d).

### Day 2 — C.1 / C.2 viral fill + E.1/E.2 map — **DONE** (Cov pull deferred)

1. [x] **C.1** paste H3N2 / COVID / HIV from existing baselines + HIV tables; Protein families + Antibody topo + H1 RF = blockers. File: `table_c1_full_tree_gen.md`. H3 TreeSBM additives filled from mrs0.5 enrich.
2. [x] **C.2** assemble Cov@100 from Table 5 / H3 eabs / HIV cov / Ab Track C; **Cov@500 = —** (blocker). File: `table_c2_future_lineage.md`.
3. [x] Map repo **Table 8** → appendix **E.1/E.2** with cons / aa|hit / antigenic. File: `table_e1_e2_ablations.md`.
4. [~] Pull Day-1 Cov jobs: **7628260/61/62** still **RUNNING** (~43m @ 12:45 ET; CSV placeholders 0 bytes) → leave D.1 Cov cells **—**; do **not** resubmit. Re-check Day 3 morning.

### Day 3 — Antibody + HIV appendix polish — **DONE** (Cov pull deferred: SSH)

1. [x] Ab Neutral SHM / pLM / AR: **REDO 2026-08-15** — Neutral SHM = **JC69** (≠ Thrifty); AR Rod.82 **FILLED**; CoSiNE/TreeSBM **FILLED**; pLM Rod.82 **in progress** (nohup). See `table6_ab_track_c.md`.
2. [x] Antibody **C.1** RF/Quartet/W1: leave **—** · **Blocker:** forced topo / no free-gen Newick (optional free-topo job **not** started).
3. [x] HIV C.1/C.2 already pasted Day 2 — captions OK; OAS **test-set** Track A row = **—** · **Blocker:** no `eval_oas_test*` artifact (Rod.82 path used).
4. [~] **Pull Cov 7628260–62:** Betty SSH from agent **timed out** (~13:00 ET; ControlMaster absent / port 22 unreachable). Leave D.1 Cov cells **—**; **do not invent**; **do not resubmit**. Re-check once Betty reachable (Day 4 morning).

### Day 4 — D.2 map + panel Table 8 (main) if required — **DONE** (Cov pull blocked)

1. [~] **Pull Cov 7628260–62:** Betty SSH **port 22 timeout** again (~13:27 ET Day 4; ping OK, `nc` fail; user interactive ssh also exit 255). Leave D.1 Cov cells **—**; **do not invent**; **do not resubmit** (status still UNKNOWN — need sacct before any resubmit).
2. [x] **D.2** mapped from NeutralBD + T7 ESM2±fit + T8 no_bridge / no_seq_branch / full → [`table_d2_ref_components.md`](table_d2_ref_components.md) + PASTE §3g. No new trains.
3. [x] Main-text **Table 8** panels: left **—** · **Blocker:** panel selector not run (author optional; deferred).
4. [x] C.3 Tree-KL: not started (optional / expensive).
5. [x] No Day 5 bulk overnight trains queued.

### Day 5 — Deferred / figures / cost / sensitivity — **DONE** (Cov pull still blocked)

1. [~] **Pull Cov 7628260–62:** Betty SSH **port 22 timeout** again (~13:31 ET Day 5; ConnectTimeout=15). Leave D.1 Cov cells **—**; **do not invent**; **do not resubmit** (status still UNKNOWN — need sacct before any resubmit). No md sync to Betty.
2. [x] **E.3** mapped from existing H3/COVID K-curves → [`table_e3_sampling.md`](table_e3_sampling.md). Temp / max-node / K=500 / sampling-horizon = **—** (no new GPU grid; Betty down).
3. [x] **F.1 / F.2 / G.1 / Fig G.1:** left **—** (deferred; no surplus compute / harness).
4. [x] **I.1** → [`table_i1_compute.tex`](table_i1_compute.tex) + [`table_i1_compute.md`](table_i1_compute.md): train **8.0 h** (covid_v5_mutrec ep14), Params **11M** reconciled; GPU mem **—**.
5. [x] **J.2 / J.3** from code/config → [`table_j2_j3_hyperparams.md`](table_j2_j3_hyperparams.md).
6. [x] **H.1 / H.2** figure selection → [`table_h1_h2_qualitative.md`](table_h1_h2_qualitative.md) (`covid_tree_viz` / `h3n2_tree_viz`).

---

## 2. Table-by-table column checklist

### A.1 Dataset construction — **DONE**

| Column | Status | Notes |
|---|---|---|
| Root assignment | DONE | PDF text |
| Leaf sequences | DONE | |
| Internal sequences | DONE | ASR / latent |
| Branch lengths | DONE | FastTree / lineage |
| Filtering | DONE | |

### A.2 Evaluation splits — **DONE**

| Column | Status |
|---|---|
| Training data | DONE |
| Test data | DONE |
| Generalization tested | DONE |

### B.1 Tree-level metrics — **DONE** (definitions)

Tree-KL, Split-KL, RF, Weighted RF, Quartet, Branch W1, Tree depth error, Imbalance error — glossary only.

### B.2 Future-lineage metrics — **DONE** (definitions)

Coverage@K,ε; Mutation recall/precision; Jaccard; Clade recall; Median min distance; Likelihood rank.

### B.3 Sequence metrics — **DONE** (definitions)

pLM NLL; Terminal edit; Motif preservation; Fitness; Diversity; Novelty.

---

### C.1 Full tree generation results

**Paper columns:** Dataset, Method, Tree-KL↓, Split-KL↓, RF↓, Quartet↓, Branch W1↓, pLM NLL↓  
**Additive (recommended):** Cons↑, Antigenic↑, AA|hit↑  
**LaTeX paste:** [`table_c1_full_tree_gen.tex`](table_c1_full_tree_gen.tex) (2026-08-17; **HIV Env → Influenza HA H1N1**; COVID/H3/Ab unchanged). Companion [`table_c1_full_tree_gen.md`](table_c1_full_tree_gen.md).  
**Prior gap jobs:** [`table_c1_job_ids.json`](table_c1_job_ids.json) **2026-08-16**. **H1N1 jobs:** [`table_c1_h1n1_job_ids.json`](table_c1_h1n1_job_ids.json) **2026-08-17** (`scripts/betty_submit_c1_h1n1.sh`).

**Best TreeSBM per dataset (pasted):**
| Dataset | Ckpt / job | Why “best” |
|---|---|---|
| SARS-CoV-2 Spike | `covid_v5_mutrec` · topo **7575017** · NLL T8 full | Only COVID TreeSBM C.1 suite; Tree-KL saturated → RF/Q/W1 from same run |
| Influenza HA H3N2 | `h3n2_v2` · `results_baselines.csv` N=16 means | Authoritative C.1 topo CSV (TABLE_METADATA); NLL **7635752** COMPLETED (`h3n2_v3_lit_hotspot` **0.463**) |
| Influenza HA H1N1 | `h1n1_v2_lit_hotspot` · geo `data/h1n1` · topo **7646375** / NLL **7646376/77** PENDING | Table 5 + C.3 locked this ckpt (ep74). **Not** `h1n1_leafholdout` / `h1n1_temporal`. All C.1 cells `--` until jobs COMPLETED |
| HIV Env | *(removed from C.1 `.tex`)* geo · bl **7626610** · NLL **7635751** | Numbers kept in `hiv_geo/` only; replaced by H1N1 |
| Antibody | OAS `ab_oas_1m_v1` · free-topo **7635753** DONE (+ NLL **7635754** / **7635755** TreeSBM NLL **0.451**) | Track C forced topo unused for C.1; Neutral CTMC row = NeutralBD free-topo |

| Dataset × Method | Tree-KL | Split-KL | RF | Quartet | W1 | pLM NLL | Cons | Antigenic | AA\|hit | Blocker / artifact |
|---|---|---|---|---|---|---|---|---|---|---|
| Protein families × all | — | — | — | — | — | — | — | — | — | **No PFAM dataset / jobs in repo** (not in `.tex`) |
| SARS-CoV-2 × Neutral | **0.693** | **43.0** | **0.980** | **0.672** | **2.29e-5** | **0.475** | — | — | — | topo `results_baselines_covid.csv` **7575017**; NLL `eval_table5_covid_baselines_plmnll.json` **7635749** |
| SARS-CoV-2 × AR | **0.693** | **62.9** | **0.990** | **0.668** | **5.62e-5** | **0.461** | **0.993** | **0.017** | **0.094** | same artifact **7635749** |
| SARS-CoV-2 × TreeSBM | **0.693** | **57.5** | **0.982** | **0.669** | **5.45e-5** | **0.457** | **0.839** | **0.005** | **0.375** | topo 7575017; NLL from T8 full / `covid_v5_mutrec` |
| Flu H3N2 × Neutral | **0.693** | **43.0** | **0.980** | **0.714** | **7.77e-4** | **0.466** | — | — | — | `results_baselines.csv` emp N=16; NLL `eval_table5_h3n2_baselines_plmnll.json` **7635750** |
| Flu H3N2 × AR | **0.693** | **62.9** | **0.986** | **0.717** | **9.21e-4** | **0.462** | — | — | — | same artifact **7635750** |
| Flu H3N2 × TreeSBM | **0.693** | **57.4** | **0.983** | **0.721** | **9.25e-4** | **0.463** | **0.918** | **0.006** | **0.439** | topo `h3n2_v2`; NLL **7635752** `…_mrs0.5_plmnll.json` |
| Flu H1N1 × Neutral | -- | -- | -- | -- | -- | -- | -- | -- | -- | topo **7646375** + NLL **7646376** PENDING; no prior `results_baselines_h1n1.csv` |
| Flu H1N1 × AR | -- | -- | -- | -- | -- | -- | -- | -- | -- | same; T5 **7539950** has AR seq but **no** `plm_nll` |
| Flu H1N1 × TreeSBM | -- | -- | -- | -- | -- | -- | -- | -- | -- | ckpt `h1n1_v2_lit_hotspot`; NLL **7646377** PENDING (mrs0.5 enrich lacks `plm_nll`) |
| Antibody × Neutral | **0.693** | **113.1** | **0.970** | **0.639** | **0.00670** | **0.513** | — | — | — | free-topo **7635753** `table_empirical_N16.csv`; NLL **7635754** |
| Antibody × AR | **0.693** | **115.4** | **0.984** | **0.648** | **0.00859** | **0.466** | — | — | — | same |
| Antibody × TreeSBM | **0.693** | **118.1** | **0.979** | **0.658** | **0.0109** | **0.451** | — | — | — | OAS ckpt; NLL **7635755** |

† Empirical Tree-KL saturates ≈ ln2 — pasted as **0.693** when a sim_neutral run exists (not left `--`). Prefer RF/Quartet/W1.  
‡ Split-KL from sim_neutral N=16 where empirical Split-KL is NaN.  
**LaTeX 2026-08-17:** HIV Env block replaced by Influenza HA H1N1 (all `--`). COVID/H3/Ab unchanged. H1 fill when **7646375–77** COMPLETED.

---

### C.2 Future-lineage coverage (multi K / ε)

**Paper columns:** Dataset, Method, Cov@100, Cov@500, Mut. recall, Clade recall, Min dist  
**Additive:** Cons, Antigenic, AA|hit

| Dataset × Method | Cov@100 | Cov@500 | Mut | Clade | Min dist | Cons / Ant / aa\|hit | Blocker |
|---|---|---|---|---|---|---|---|
| Protein families | — | — | — | — | — | — | no data |
| SARS-CoV-2 pLM / TreeSBM | **0.775 / 0.788** | pLM **0.775**; TreeSBM — | **0.048 / 0.080** | **0 / 0** | **1.425 / 1.388** | 0.985/0 / 0.017 vs 0.839/0.005/0.375 | TreeSBM Cov@500 **7635911** RUNNING (CSV 0 bytes) |
| Flu HA H3N2 pLM / TreeSBM | **0.675 / 0.725** | **0.750 / 0.788** | **0.004 / 0.100** | **0.000 / 0.100** | **2.650 / 2.238** | H3 TreeSBM enrich; pLM ant/aa — | K500 rescored e-abs |
| HIV geo TreeSBM / pLM | **0.188 / 0.188** | **0.200 / 0.188** | **0.026 / 0.015** | **0.061 / 0.064** | **63.1 / 63.9** | enrich additives | pLM **7635912–14** DONE |
| Antibody pLM / `treesbm_ab_oas` | **0.018 / 0.027** | **0.018 / 0.027** | CDR **0.297 / 0.150** | **0.039 / 0.018** | **15.9 / 14.3** | — | **7639280** eval DONE (fix comma-export on 7635917) |

---

### C.3 Forecasting by horizon — TreeSBM numbers (Tree-KL —)

**Paper columns:** Dataset, Horizon, Cov@100↑, Mut. recall↑, Clade recall↑, Min dist↓, Tree-KL↓  
**Additive:** Cons, AA|hit, Antigenic (`lit_hotspot_mut_frac`; job **7628578**)  
**Cov@100** = `coverage_obs_e2` @ K=100. Jobs **7627191–95** COMPLETED. Empty buckets locked empty (protocol). Tree-KL **—** (deferred).  
**LaTeX:** [`table_c3_horizon.tex`](table_c3_horizon.tex) (2026-08-17; user template; **only `--` replaced**). Clade recall filled from TreeSBM K=100 CSVs: COVID med/long **0.000**, H3 med/long **0.000**, HIV geo short **0.000**, HIV geo long **0.097**. Remaining `--` = empty locked buckets (COVID short, H3 short, HIV geo medium) — **not** missing jobs. HIV paper row = **geo** (matches short 0.313). Paper COVID long Mut **0.003** is NeutralBD (`0.00347`); TreeSBM long Mut is **0.000** (left as user’s already-filled cell).

| Dataset | Horizon | Method | Cov@100 | Mut | Clade | Min dist | Tree-KL | Cons | AA\|hit | Antigenic |
|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|
| SARS-CoV-2 | short | — | — | — | — | — | — | — | — | — |
| SARS-CoV-2 | medium | NeutralBD | 0.938 | 0.000 | 0.000 | 0.69 | — | 1.000 | — | 0.010 |
| SARS-CoV-2 | medium | AR | 0.938 | 0.000 | 0.000 | 0.69 | — | 1.000 | — | 0.008 |
| SARS-CoV-2 | medium | TreeSBM | 0.938 | 0.000 | 0.000 | 0.69 | — | 1.000 | — | 0.006 |
| SARS-CoV-2 | long | NeutralBD | 0.667 | 0.003 | 0.000 | 1.90 | — | 1.000 | 1.000† | 0.013 |
| SARS-CoV-2 | long | AR | 0.667 | 0.000 | 0.000 | 1.92 | — | 1.000 | — | 0.005 |
| SARS-CoV-2 | long | TreeSBM | 0.667 | 0.000 | 0.000 | 1.92 | — | 1.000 | — | 0.017 |
| Flu H1N1 | short | NeutralBD | 1.000 | 0.045 | 0.750 | 0.44 | — | 1.000 | 1.000 | 0.115 |
| Flu H1N1 | short | AR | 1.000 | 0.067 | 0.250 | 0.46 | — | 1.000 | 1.000 | 0.110 |
| Flu H1N1 | short | TreeSBM | 1.000 | 0.212 | 0.500 | 0.31 | — | 1.000 | 1.000 | 0.056 |
| Flu H1N1 | medium | — | — | — | — | — | — | — | — | — |
| Flu H1N1 | long | NeutralBD | 0.406 | 0.013 | 0.000 | 2.91 | — | 1.000 | 1.000 | 0.107 |
| Flu H1N1 | long | AR | 0.406 | 0.005 | 0.000 | 2.91 | — | 1.000 | 1.000 | 0.093 |
| Flu H1N1 | long | TreeSBM | 0.500 | 0.141 | 0.100 | 2.41 | — | 1.000 | 1.000 | 0.097 |
| Flu H3N2 | short | — | — | — | — | — | — | — | — | — |
| Flu H3N2 | medium | NeutralBD | 0.828 | 0.000 | 0.000 | 1.27 | — | 1.000 | — | 0.013 |
| Flu H3N2 | medium | AR | 0.859 | 0.014 | 0.250 | 1.22 | — | 1.000 | 1.000 | 0.026 |
| Flu H3N2 | medium | TreeSBM | 0.859 | 0.099 | 0.000 | 1.14 | — | 1.000 | 1.000 | 0.011 |
| Flu H3N2 | long | NeutralBD | 0.000 | 0.000 | 0.000 | 8.25 | — | 1.000 | — | 0.016 |
| Flu H3N2 | long | AR | 0.000 | 0.000 | 0.000 | 8.25 | — | 1.000 | — | 0.020 |
| Flu H3N2 | long | TreeSBM | 0.000 | 0.131 | 0.000 | 7.19 | — | 1.000 | 1.000 | 0.005 |
| HIV geo | short | NeutralBD | 0.313 | 0.001 | 0.067 | 88.6 | — | 1.000 | 0.523 | 0.188 |
| HIV geo | short | TreeSBM | 0.313 | 0.001 | 0.000 | 88.6 | — | 1.000 | 0.740 | 0.226 |
| HIV geo | medium | — | — | — | — | — | — | — | — | — |
| HIV geo | long | NeutralBD | 0.000 | 0.023 | 0.097 | 26.8 | — | 0.982 | 0.576 | 0.142 |
| HIV geo | long | TreeSBM | 0.000 | 0.050 | 0.097 | 24.7 | — | 0.994 | 0.695 | 0.146 |
| HIV temporal | short | — | — | — | — | — | — | — | — | — |
| HIV temporal | medium | NeutralBD | 0.000 | 0.016 | 0.000 | 177.6 | — | 1.000 | 0.694 | 0.212 |
| HIV temporal | medium | TreeSBM | 0.000 | 0.020 | 0.113 | 177.1 | — | 1.000 | 0.701 | 0.191 |
| HIV temporal | long | NeutralBD | 0.000 | 0.023 | 0.067 | 186.5 | — | 0.973 | 0.212 | 0.142 |
| HIV temporal | long | TreeSBM | 0.000 | 0.037 | 0.098 | 184.3 | — | 0.994 | 0.400 | 0.135 |

† sparse / nan-prone when almost no generated mutations hit GT sites. Empty buckets: COVID short; H1 med; H3 short; HIV geo med; HIV temporal short (locked T5/T4 roots).  
Paste file: `benchmarks/results/tables/table_c3_horizon.md`.

---

### D.1 Full pLM reference-process comparison

**Paper columns:** Reference, Params, Fitness, Tree-KL↓, Split-KL↓, Cov@100↑, Mut. recall↑, pLM NLL↓  
**Additive:** Cons↑, Antigenic↑, AA|hit↑

| Reference | Params | Fitness | Tree-KL | Split-KL | Cov@100 | Mut | pLM NLL | Cons | Ant | AA\|hit | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Substitution-only | — | no | — | — | **0.775** | **0.022** | **0.471** | **0.821** | **0.007** | **0.083** | DONE (JTT) |
| ESM-2-650M | 650M | no | — | — | **0.775** | **0.018** | **0.442** | **0.893** | **0.006** | **0.108** | DONE |
| ESM-2-650M | 650M | yes | — | — | **0.775** | **0.001** | **0.441** | **0.988** | **0.002** | **0.098** | DONE |
| ESM-C | — | no | — | — | **0.775** | **0.035** | **0.461** | **0.692** | **0.006** | **0.080** | Cov **7628260** DONE |
| ESM-C | — | yes | — | — | **0.775** | **0.036** | **0.427** | **0.713** | **0.005** | **0.085** | Cov **7628261** DONE |
| ProGen2 | — | yes | — | — | **0.788** | **0.113** | **0.455** | **0.985** | **0.005** | **0.398** | Cov **7628262** DONE (e2=0.7875) |
| Evo2 | — | yes | — | — | — | — | — | — | — | — | **Blocker: Evo2 backend not wired** |
| TreeSBM default | — | yes | 0.693† | 57.52 | **0.775** | **0.038** | **0.446** | **0.941** | **0.006** | **0.290** | DONE |

Per-R0 Tree-KL/Split-KL: **blocker** (would need `MODE=baselines` per ROW; Tree-KL still saturates).

---

### D.2 Reference process component ablation — **MAPPED Day 4**

**Paper columns:** Variant flags + Cov@100↑, Tree-KL↓, Branch W1↓  
**Additive:** Cons, Antigenic, AA|hit (+ Mut / pLM NLL)  
**Paste:** [`table_d2_ref_components.md`](table_d2_ref_components.md)

| Variant | Map | Cov@100 | Tree-KL | Branch W1 | Cons | Ant | AA\|hit |
|---|---|---:|---:|---:|---:|---:|---:|
| Neutral birth-death | NeutralBD T5 eabs + **7575017** | 0.775 | 0.693† | 2.29e-5 | 1.000 | **0.012** | **0.828** |
| pLM mutation only | T7 ESM2 nofit | 0.775 | — | **5.20e-5** | 0.893 | 0.006 | 0.108 |
| pLM + fitness | T7 ESM2 fit | 0.775 | — | **5.32e-5** | 0.988 | 0.002 | 0.098 |
| pLM + branching | T8 full (β=0; vs `no_seq_branch` OFF) | 0.775 | 0.693† | 5.26e-5 | 0.839 | 0.005 | 0.375 |
| pLM + fitness + branching | T8 no_bridge (**careful caption**) | 0.775 | 0.693† | 5.49e-5 | 0.876 | 0.004 | 0.086 |
| TreeSBM full bridge | T8 full | 0.775 | 0.693† | 5.26e-5 | 0.839 | 0.005 | 0.375 |

**Blockers cleared:** ESM2 Branch W1 **7628576/77** COMPLETED. Remaining: dedicated factorial not run; ESM2 Tree-KL — (saturated); main Table 8 panels —.

---

### E.1 Bridge matching ablations — **PARTIAL 2026-08-17**

**Paper columns:** Bridge matching, Doob, Terminal loss flags + Tree-KL↓, Cov@100↑, Mut. recall↑  
**Additive:** Cons, Antigenic, AA|hit, pLM NLL  
**Paste:** [`table_e1_e2_ablations.md`](table_e1_e2_ablations.md) · jobs [`table_e1_job_ids.json`](table_e1_job_ids.json)  
**Protocol:** COVID Brazil / `covid_v5_mutrec` recipe (same as T8). Tree-KL **—** on new rows (would saturate ln2; not queued).

| Variant | Map | Tree-KL | Cov@100 | Mut | Cons | Ant | AA\|hit | pLM NLL |
|---|---|---|---:|---:|---:|---:|---:|---:|
| Reference process only | T8 no_bridge / T7 R0 | 0.693 | 0.775 | 0.008 | 0.876 | 0.004 | 0.086 | 0.4465 |
| Terminal-only generator | `--ablate-terminal-only` **7628589** / en **7628592** / cov **7628593** | — | **0.775** | **0.089** | **0.797** | **0.009** | **0.231** | **0.4609** |
| Bridge w/o Doob form | `--ablate-doob` **7628590** (+en **7628594** / cov **7628595**) | — | — | — | — | — | — | — |
| TreeSBM w/o terminal consistency | `--ablate-terminal-consistency` **7628591** / en **7628596** / cov **7628597** | — | **0.775** | **0.062** | **0.794** | **0.008** | **0.201** | **0.4599** |
| TreeSBM full | T8 full | 0.693 | 0.775 | 0.080 | 0.839 | 0.005 | 0.375 | 0.4572 |

**Status:** trains **7628589** / **7628591** COMPLETED + enrich/cov COMPLETED (pulled). **7628590** no_doob still **RUNNING** (~1-21h of 72h); afterok **7628594/95** PENDING Dependency. Do not invent no_doob metrics.

---

### E.2 Architecture ablations — **QUEUED 2026-08-17**

**Paper columns:** Tree encoder / Branch / Mutation / Stop heads + Tree-KL↓, Cov@100↑  
**Additive:** Cons, Antigenic, AA|hit, Mut, pLM NLL (Cov@100 + pLM NLL kept on top of Tree-KL)

| Variant | Map / flag | Cov@100 | Tree-KL | Mut | Cons | Ant | AA\|hit | pLM NLL |
|---|---|---:|---|---:|---:|---:|---:|---:|
| Node-independent rates | T8 no_tree_ctx | 0.775 | — | 0.084 | 0.845 | 0.012 | 0.264 | 0.4562 |
| No branch-length head | T8 no_bl_eval | 0.775 | 0.693 | 0.080 | 0.839 | 0.005 | 0.375 | 0.4572 |
| No mutation head | `--ablate-mut-head` train **7642530** (+en **7642532** / cov **7642533**) | — | — | — | — | — | — | — |
| No stop/termination head | `--ablate-stop-head` train **7642531** (+en **7642534** / cov **7642535**) | — | — | — | — | — | — | — |
| Full TreeSBM | T8 full | 0.775 | 0.693 | 0.080 | 0.839 | 0.005 | 0.375 | 0.4572 |

72h train wall; `--resume` in `checkpoints/covid_e2_{no_mut_head,no_stop_head}/`. IDs: [`table_e2_job_ids.json`](table_e2_job_ids.json). Metric cells stay **—** until afterok enrich/cov land.

λ_br=0 retrain **7539981** TIMEOUT — do not use as finished row.

---

### E.3 Sampling hyperparameter sensitivity — **TEMP/MAXL FILLED 2026-08-17**

**Paste:** [`table_e3_sampling.md`](table_e3_sampling.md) · IDs [`table_e3_job_ids.json`](table_e3_job_ids.json)

| Setting | Value | Coverage | Diversity | Tree-KL | Runtime | Status |
|---|---|---|---|---|---|---|
| K=50/100 (H3 e5 / COVID e2) | 50, 100 | H3 e5 **0.812/0.812**; COVID e2 **0.775/0.788** | — | — | 15668 / 21421 | Existing curves |
| K=500 (COVID e2) | 500 | — | — | — | — | **7635911** RUNNING (CSV 0 bytes) |
| Temperature | 0.5 / 1.0 / 1.5 | **0.775** / **0.788** / **0.775** | — | — | 20498 / 21421 / 20639 | **7639276/77** COMPLETED; 1.0=Table5 |
| Max leaves | 100 / 400 / 800 | **0.775** / **0.788** / **0.775** | — | — | 24655 / 21421 / 44789 | **7639278/79** COMPLETED; 400=prod default |
| Evolutionary horizon | — | — | — | — | — | C.3 = genetic-H, not sampling-horizon grid |

---

### F.1 Data scaling — **MISSING / DEFER**

| Training fraction | # roots | # trees | Tree-KL | Cov@1000 | Mut | pLM NLL |
|---|---|---|---|---|---|---|
| 10/25/50/75/100% | — | — | — | — | — | — |

**Blocker:** no fractioned retrain schedule; Cov@1000 Track A unfinished.

---

### F.2 Cross-domain transfer — **MISSING / DEFER**

All train→test domain cells **—**.  
**Blocker:** no cross-domain eval harness locked to paper rows.

---

### G.1 Likelihood ranking — **MISSING / DEFER**

| Dataset | True median rank | Random decoy | Matched-mutation decoy | Shuffled-path decoy |
|---|---|---|---|---|
| All rows | — | — | — | — |

**Blocker:** ranking/decoy script not producing paper table artifacts.

### Fig G.1 Coverage calibration — **MISSING / DEFER**

---

### H.1 / H.2 Qualitative — **DOCUMENTED Day 5**

**Paste:** [`table_h1_h2_qualitative.md`](table_h1_h2_qualitative.md). COVID g40 matched + H3 g55/g40 screen ranks; HIV/Ab galleries **—**.

---

### I.1 Computational cost — **FILLED 2026-08-17**

**Paste:** [`table_i1_compute.tex`](table_i1_compute.tex) · notes [`table_i1_compute.md`](table_i1_compute.md)

| Method | Params | Train time | GPU mem | Time/root | Notes |
|---|---|---|---|---|---|
| pLM prior | ESM-2 8M | --- (no train) | — | 3 s COVID | Table 5 coverage runtime |
| Autoregressive | — | — | — | 31 s | AR params + train wall — (external) |
| TreeSBM | **11M** (~3M train + 8M PLM) | **8.0 h** to ep14 (`covid_v5_mutrec`, job **7350224**) | — | 71 min COVID | `runtime_gen_sec`=21421 / 5 roots |

**Reconciled Params:** paper **11M** = ~2.9M trainable in `best.pt` + 8M frozen ESM-2 (not the old “~3M only” shorthand). **GPU mem** still — (no sacct MaxRSS / nvidia-smi pull).

---

### J.1 Linear-probe validation — **DONE (FINAL)**

**Paper columns:** Embedding × {Leaf/Int AUROC, P–C AUROC, AncAA F1, MutPos F1, Subtree R², Depth R², RootDist R², Subs R², HamPar R²}

| Embedding | Status | Artifact |
|---|---|---|
| ESM (raw) | DONE | wave4 `transformer_val_*_j1` |
| ESM+branch | DONE | |
| Topology | DONE | |
| Graphrand | DONE | |
| NodeEnc | DONE | |
| GraphTF | DONE | |

Primary paste already in `PAPER_TABLES.md` (H3N2 + COVID). Paper table is H3N2-shaped; COVID companion optional.

---

### J.2 Model hyperparameters — **DONE Day 5; LaTeX filled 2026-08-17**

**Paste:** [`table_j2_j3_hyperparams.tex`](table_j2_j3_hyperparams.tex) · notes [`table_j2_j3_hyperparams.md`](table_j2_j3_hyperparams.md)

| Hyperparameter | Value | Source |
|---|---|---|
| Tree encoder layers | **4** | `train.py` |
| Hidden dimension | **128** | `train.py` |
| Attention heads | **8** | `train.py` |
| pLM reference | **ESM-2 8M** | `esm2_t6_8M_UR50D` |
| Mutation / branching / BL head depth | **2 / 2 / 2** | `networks.py` (`deep_mut_head=false`) |
| Dropout | **0.1** | `train.py` |
| Optimizer / LR / batch | **AdamW / 1e-4 / 1** | `train.py` + slurm |
| Training epochs | **14** (best ckpt) | job **7350224** / `covid_v5_mutrec/best.pt` |

### J.3 Sampling settings — **DONE Day 5; LaTeX filled 2026-08-17**

Same files. Table 5 Coverage@100 protocol (`covid_v5_mutrec`, job **7575015**):

| Setting | Value | Notes |
|---|---|---|
| Evolutionary horizon | **50** | `n_steps`; enrichment uses **100** |
| Maximum nodes | **400** | `max_leaves` prod default |
| Maximum branching factor | **2** | binary bifurcation |
| Mutation temperature | **1.0** | E.3 flat |
| Branching temperature | **1.0** (learned $\lambda$; scale **6.0**) | no temp knob; `branch_rate_scale=6.0` |
| Stop threshold | **none** (fixed `n_steps`) | stop head not used at inference |
| K | **100** | Cov@100 |
| Seeds | **0** | coverage; enrich **42** |

---

### K Algorithms — **DONE** (pseudo-code in PDF)

---

## Coverage flatness (2026-08-15)

On locked COVID Brazil 5-roots, **Cov@ε=2 @ K=100 ≈ 0.775** for NeutralBD / pLM±fit / all Table-8 ablations: easy leaves saturate and one hard root stays uncovered. **ε, virus, and horizon** move coverage; bridge/entropy/mask mainly move Mut / aa|hit. Full writeup: [`COVERAGE_WHAT_MOVES_IT.md`](COVERAGE_WHAT_MOVES_IT.md).


| Item | Result |
|---|---|
| Job **7627915** (OAS precompute) | **COMPLETED** 00:04:03 — all **complete** groups cached |
| Complete groups with `*_plm.pt` + `*_ref_rates.pt` | train **397**, val **44**, test **170** (matches donor-disjoint clone set) |
| Extra `group_*.nwk` without `*_rooted.nwk`+`*_anc_aa.fasta` | **not** precompute targets (`precompute_plm.py` only indexes complete groups) |
| Job **7627919** (Rod.82 eval) | **COMPLETED** 01:03:40 → `treesbm_ab_oas` metrics in `table6_ab_track_c.*` |
| **Is OAS PLM precompute still required?** | **No** for Rod.82 eval (uses ckpt + observed Rod.82 trees). Precompute **already finished** for OAS complete groups. Only needed again if adding new incomplete clones into training. |

OAS paste row (Track C): Cov@e2=**0.024**, SHM err=**0.077**, CDR recall=**0.150**, term-div err=**9.89**, RF=NaN (forced topo).

### Track A v2 — SHM Q0 retrain (submitted 2026-08-20)

v1 CDR=0.150 is a **Q0 mismatch** (ESM-2 MLM vs SHM targeting), not a missing viral mut-recovery flag. v2 injects CDR ∪ AID/WRCH stay-logit boost into Q0 (`--shm-site-boost 2.0 --shm-fwr-stay 0.5 --shm-use-aid`); ckpt `checkpoints/ab_oas_1m_v2_shm`; Rod samples `treesbm_ab_oas_v2`. Does not overwrite v1 ckpt/samples or pathogen `best.pt`. Same 397/44/170 splits; PLM from **7627915**. Jobs: train **7720919**, Rod.82 **7720920** — [`ab_oas_v2_shm_job_ids.json`](ab_oas_v2_shm_job_ids.json). Entropy audit **7720923** COMPLETED (OAS H_CDR/H_FWR=1.627, Rod=2.024). Metrics **—** until Rod eval completes.

Full wave tracker (viral eval/retrain, VaxSeer, flu path): [`RETRAIN_VAXSEER_SHM_WAVE.md`](RETRAIN_VAXSEER_SHM_WAVE.md).


---

## 4. Horizon C.3 (B) — locked answer

| Job | Virus | State |
|---:|---|---|
| 7627191 | covid | COMPLETED |
| 7627192 | h1n1 | COMPLETED |
| 7627193 | h3n2 | COMPLETED |
| 7627194 | hiv_geo | COMPLETED |
| 7627195 | hiv_temporal | COMPLETED |

**Next wave started:** local aggregation → `table_c3_horizon.md` + job JSON update.  
**Not started (protocol forbids / expensive):** empty-bucket root expansion; pLM-prior C.3 methods; Tree-KL.

---

## 5. Jobs started this session vs deferred

### Started / Day 1 completed (2026-08-15 ~12:00 ET)
- [x] C.3 paste (`table_c3_horizon.md` + `PASTE_TABLES` §3c)
- [x] Ab OAS Rod.82 `treesbm_ab_oas` paste (`table6_ab_track_c.md` + `PASTE_TABLES` §3b)
- [x] ProGen2 enrich **7627918** n=20 → `table7_refproc.md` / PASTE §2 (cons / aa|hit / antigenic / mut / NLL)
- [x] J.1 FINAL confirmed (`PAPER_TABLES.md` / PASTE §3d)
- T7 Cov@100 overnight (still running — **do not resubmit**): **7628260** / **7628261** / **7628262**

### Day 2 completed (2026-08-15 ~12:45 ET)
- [x] C.1 viral paste + additives → `table_c1_full_tree_gen.md` + PASTE §3 (H3 TreeSBM cons/ant/aa|hit from mrs0.5 enrich)
- [x] C.2 assemble → `table_c2_future_lineage.md` + PASTE §3e (Cov@500 —)
- [x] E.1/E.2 map from T8 → `table_e1_e2_ablations.md` + PASTE §3f
- [x] Checked Cov **7628260–62**: still **RUNNING** (~43m); leave D.1 Cov — ; no new overnight trains started

### Day 2 remaining / carry to Day 3
- Cov@100 for ESM-C ±fit + ProGen2 still **RUNNING** (7628260–62) → pull when COMPLETED
- Evo2 D.1 row still **—** (backend not wired)
- C.3 Tree-KL / antigenic on C.3 gens still **—**
- Ab Neutral SHM / pLM / AR Rod.82 caption; Ab C.1 free-topo RF optional
- H3 Neutral/AR seq enrich; Protein-family C.1/C.2; Cov@500; **H1N1 C.1 queued 7646375–77**

### Day 3 completed (2026-08-15 ~13:00 ET) + Table 5 redo (~14:00 ET)
- [x] Ab Neutral SHM caption map (thrifty ≡ Neutral SHM) + Ab pLM/AR — blockers → `table6_ab_track_c.md` / `antibody_tracks.md` / PASTE §3b *(superseded by redo below)*
- [x] **Table 5 redo:** Neutral SHM = **JC69** Rod.82 (≠ Thrifty); AR Rod.82 **FILLED**; CoSiNE CDR/SHM **FILLED**; TreeSBM pathogen **FILLED**; pLM Rod.82 **RUNNING** (nohup; sbatch DNS broken on login01)
- [x] Ab C.1 free-topo RF left **—** (no job); OAS Track A test-set left **—** (no artifact)
- [~] Cov **7628260–62**: Betty SSH **timeout** — could not confirm COMPLETED/FAILED/RUNNING; leave D.1 Cov **—**; no resubmit; no Day 4–5 bulk trains started

### Day 4 completed (2026-08-15 ~13:30 ET)
- [x] **D.2** map → `table_d2_ref_components.md` + PASTE §3g (NeutralBD / ESM2±fit / T8 no_bridge / full; additives on top)
- [x] Main Table 8 panels left **—** (panel selector not run)
- [~] Cov **7628260–62**: Betty SSH **port 22 timeout** — still cannot sacct/scp; leave D.1 Cov **—**; no resubmit; no Day 5 bulk trains started

### Day 5 completed (2026-08-15 ~13:35 ET)
- [x] **E.3** → `table_e3_sampling.md` + PASTE §3h (K=50/100 from existing curves; temp/K500/max-node —)
- [x] **I.1** → `table_i1_compute.tex` + `table_i1_compute.md` + PASTE §3i (train 8.0 h; GPU mem —)
- [x] **J.2 / J.3** → `table_j2_j3_hyperparams.md` + PASTE §3j
- [x] **H.1 / H.2** → `table_h1_h2_qualitative.md` + PASTE §3k
- [x] F.1 / F.2 / G.1 / Fig G.1 left **—** (deferred; no new trains)
- [~] Cov **7628260–62**: Betty SSH **port 22 timeout** — still cannot sacct/scp; leave D.1 Cov **—**; no resubmit; **no Betty md sync**

### Remaining open blockers (post–Day 5 + coverage session)
- **C.1 H1N1** (HIV Env replaced): topo **7646375** + Neutral/AR NLL **7646376** + TreeSBM NLL **7646377** PENDING — leave `--`; do **not** invent; pull `results_baselines_h1n1.csv` / `eval_table5_h1n1_baselines_plmnll.json` / `eval_enrichment_h1n1_v2_lit_hotspot_mrs0.5_plmnll.json`
- **Cov@500 C.2** gap jobs queued 2026-08-16 login03 — do **not** invent; pull when COMPLETED: COVID TSBM **7635911**; HIV pLM **7635912–14**; Ab **7635915–17** (see `table_c2_cov500_job_ids.json`)
- **C.3 antigenic** score job **7628578** — **FILLED** in `table_c3_horizon.md` / APPENDIX C.3
- **D.2** ESM2 Branch W1 **7628576/77**; Neutral ant/aa from **7628578**
- Evo2 D.1 row **—** (backend not wired)
- C.3 Tree-KL **—**
- E.1 terminal-only **FILLED** (Cov **0.775**, Mut **0.089**, Cons **0.797**, Ant **0.009**, aa|hit **0.231**, NLL **0.4609**); w/o terminal **FILLED** (Cov **0.775**, Mut **0.062**, Cons **0.794**, Ant **0.008**, aa|hit **0.201**, NLL **0.4599**); w/o Doob train **7628590** still RUNNING + afterok **7628594/95** PENDING
- E.2 no mut/stop heads **queued** train **7642530/31** (72h) + afterok enrich/cov **7642532–35**; metrics —
- E.3 temp/maxl **FILLED** (**7639276–79**); K=500 **7635911** RUNNING (CSV 0 bytes); Diversity —
- F/G entirely
- Coverage flatness writeup: [`COVERAGE_WHAT_MOVES_IT.md`](COVERAGE_WHAT_MOVES_IT.md)
- **Ab OAS v2 SHM retrain** submitted 2026-08-20: train **7720919** + afterok Rod.82 **7720920** → `treesbm_ab_oas_v2` (see `ab_oas_v2_shm_job_ids.json`). Metrics — until COMPLETED. Do not overwrite v1 `treesbm_ab_oas` / `ab_oas_1m_v1`.
- **Viral mutlin wave** eval **7720891–910** + trains **7720911–918** (β∈{0,0.25}); tracker [`RETRAIN_VAXSEER_SHM_WAVE.md`](RETRAIN_VAXSEER_SHM_WAVE.md).

**Coverage note (paste):** On locked COVID Brazil 5-roots, Cov@ε=2@K=100 stays ≈0.775 across NeutralBD / pLM±fit / all T8 ablations because easy leaves saturate and one hard root stays uncovered; ε, virus, and horizon move coverage — bridge/entropy/mask mostly move Mut/aa|hit. See `COVERAGE_WHAT_MOVES_IT.md`.

### Deferred (explicitly)
- C.3 Tree-KL; empty-bucket fills; C.3 antigenic enrichment
- Evo2 D.1 row; F.1/F.2; G.1; E.3 full factorial grid; main Table 8 panels
- Ab Neutral SHM / pLM / AR as dedicated Rod.82 generators (**JC69 Neutral + AR DONE 2026-08-15**; **pLM Rod.82 still running** / sbatch via login03 when DNS works)
- Protein-family C.1/C.2 entirely
- New COVID E.3 GPU grid (blocked by Betty + Day 5 scope)
---

## 6. Additive metric policy (must not drop existing columns)

Wherever viral/enrichment tables are pasted (C.1 seq cols, C.2, C.3, D.1, E.*, main T5/T7/T8):

1. Keep paper’s original columns.
2. **Add** (when available in JSON/CSV):
   - `cons_retention`
   - antigenic site proportion (`lit_hotspot_mut_frac` / flu lit / HIV V1–V5)
   - `aa_acc_given_hit` (aa|hit)
3. If additive missing: write **—** + blocker (e.g. “C.3 runner lacks antigenic”).

Leaf vs tree scoring: see `EVAL_METRICS_LEAF_VS_TREE.md` (primaries are leaf-best-match).

### COVID pLM EVEscape Δ gap (main-text Table 4 / repo Table 5)

| Cell | Status | Blocker / fix |
|---|---|---|
| SARS-CoV-2 × pLM prior × EVEscape Δ | **—** (do not invent) | **7539949** `eval_table5_covid_baselines.json`: `plm_prior` EVEscape = **NaN** — ESM L=566 OOB on COVID L=1280 (17/20 trees errored). AR (−0.371) and TreeSBM (−0.048) OK from same / T8 enrich. |
| Re-run | **7644883** PENDING | `sbatch --export=ALL,VIRUS=covid,METHODS=plm_prior,OUT=checkpoints/eval_table5_covid_plm_evescape.json scripts/slurm_table5_baselines.sh` — pull `evescape_score_delta` when COMPLETED → fill `table5_viral_forecasting.{md,csv}` + PASTE §3l. |
