# Appendix D.2 — Reference-process component ablation (mapped Day 4)

**Updated:** 2026-08-15 Day 4.  
**Domain:** COVID Brazil geo-test (same roots / ckpt as D.1 / E.*).  
**Rule:** map overlapping T7 / T8 / NeutralBD cells — **no new trains**. Caption as mapped factorial, not a dedicated design. Empty cells stay **—** with blocker.

**Paper columns:** Cov@100↑, Tree-KL↓, Branch W1↓  
**Additive (on top):** Cons↑, Antigenic↑ (PMC lit), AA|hit↑ (+ Mut / pLM NLL where available)

Cov@100 = `coverage_obs_e2` @ K=100. Same 5 Brazil roots → **≈0.775** for almost every method (weak discriminator). Prefer Mut / aa|hit / min_edit / Cons.

Tree-KL† = sim_neutral ≈ ln2 (saturated). Branch W1 from baselines / T8b when present; **—** for per-R0 T7 backends (not computed).

---

## Paste-ready

| Variant | Map (artifact) | Cov@100 | Tree-KL | Branch W1 | Cons | Antigenic | AA\|hit | Mut | pLM NLL |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Neutral birth-death | NeutralBD `coverage_curves_covid_N16_table5_eabs.csv` K=100 + baselines **7575017**; ant/aa|hit cache score **7628578** → `table_d2_neutralbd_antigenic.json` | **0.775** | 0.693† | **2.29e-5**‡ | **1.000**§ | **0.012**¶ | **0.828**¶ | 0.000§ / 0.005¶ | — |
| pLM mutation only | T7 ESM-2-650M nofit (`--ablate-bridge`, β=0) enrich **7626623** / cov **7585985**; Branch W1 job **7628576** | **0.775** | — | **5.20e-5** | **0.893** | **0.006** | **0.108** | 0.018 | **0.442** |
| pLM + fitness | T7 ESM-2-650M fit (`--ablate-bridge`, β=1) enrich **7626624** / cov **7585986**; Branch W1 job **7628577** | **0.775** | — | **5.32e-5** | **0.988** | **0.002** | **0.098** | 0.001 | **0.441** |
| pLM + branching | T8 **full** (β unset=0; learned branching + bridge) vs OFF companion `no_seq_branch` | **0.775** | 0.693† | **5.26e-5** | **0.839** | **0.005** | **0.375** | 0.080 | **0.457** |
| pLM + fitness + branching | T8 `no_bridge` (ablate-bridge) — **careful caption** | **0.775** | 0.693† | **5.49e-5** | **0.876** | **0.004** | **0.086** | 0.008 | **0.447** |
| TreeSBM full bridge | T8 **full** (same row as “pLM + branching” for Cov/topo; seq enrich matches) | **0.775** | 0.693† | **5.26e-5** | **0.839** | **0.005** | **0.375** | 0.080 | **0.457** |

† Tree-KL saturated (sim_neutral ln2). ‡ Empirical N=16 mean `branch_w_all` for NeutralBD in `results_baselines_covid.csv` job **7575017** (matches C.1). § From coverage CSV (cons=1.0; mut≈0). ¶ Cache-scored from Table5 NeutralBD trees (**7628578**): ant=`lit_hotspot_mut_frac` **0.0115**; aa|hit **0.828** (2/5 roots had finite aa|hit; mean of those). Branch W1 ESM2±fit from **7628576/77** COMPLETED → `results_baselines_covid_d2_esm2_{nofit,fit}.csv` (N=16 mean `branch_w_all`).

### New jobs (2026-08-15 → pull 2026-08-16)

| Piece | Job | Status |
|---|---:|---|
| ESM2 nofit Branch W1 | **7628576** | **COMPLETED** → W1 **5.20e-5** |
| ESM2 fit Branch W1 | **7628577** | **COMPLETED** → W1 **5.32e-5** |
| NeutralBD ant + aa\|hit + C.3 ant | **7628578** | COMPLETED 00:03:45 |

See `table_d2_job_ids.json`.

### Branching OFF companion (not a paper D.2 row; caption contrast)

| T8 ablation | Cov@100 | Tree-KL | Branch W1 | Cons | Ant | AA\|hit | Mut | pLM NLL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w/o seq-dependent branching (`no_seq_branch`, jobs enrich **7597321** / cov **7597327**) | 0.775 | —¶ | —¶ | 0.818 | 0.011 | 0.242 | 0.038 | 0.454 |

¶ **Blocker:** T8b KL/W1 jobs not run for `no_seq_branch` (would be ln2 / near-full W1 anyway).

---

## Mapping notes / honesty

1. **Not a true factorial.** Overlapping T7 R0 backends (ablate-bridge) and T8 ablations; caption *“mapped from D.1 / Table 8 / NeutralBD baselines.”*
2. **pLM + branching** uses T8 **full** (learned branching ON, fitness β=0). Isolated “branching-only without bridge” does **not** exist. Contrast OFF state = `no_seq_branch`.
3. **pLM + fitness + branching** maps to T8 **`no_bridge`** per plan (“treesbm-ish without bridge”) with **careful caption:** generation is R0 / ablate-bridge (mut/aa|hit look like reference-only, not full TreeSBM heads). Do **not** claim a learned-branching + fitness − bridge factorial.
4. **TreeSBM full bridge** = T8 full (same Cov/Tree-KL/W1/additives as “pLM + branching” row). For fitness-on TreeSBM, see D.1 TreeSBM β=1 (Cov 0.775, Cons 0.941, aa|hit 0.290, Mut 0.038) — optional caption cross-ref, not a separate D.2 paper column.
5. **Per-R0 Tree-KL** for ESM2 rows = **—** (skip; saturates). **Branch W1** jobs **7628576/77** COMPLETED — ESM2 nofit **5.20e-5**, fit **5.32e-5**.
6. **Main-text Table 8 panels** (Diversity / Future cov / Antigenic cov / Worst-case / Redundancy) = **—** · **Blocker:** panel selector not run (Day 4 optional; deferred).

---

## Job / artifact cheat sheet

| Piece | Job / file |
|---|---|
| NeutralBD Cov@100 | `coverage_curves_covid_N16_table5_eabs.csv` (K=100) |
| NeutralBD Tree-KL / W1 | `results_baselines_covid.csv` **7575017** |
| ESM2 ±fit | enrich **7626623/24**; cov **7585985/86**; `table7_refproc.md` |
| T8 no_bridge | enrich **7597319**; cov **7585972**; KL **7585977** → `results_baselines_covid_t8_no_bridge.csv` |
| T8 no_seq_branch | enrich **7597321**; cov **7597327**; KL **—** |
| T8 full | enrich **7585970**; cov **7585971**; KL **7585976** → `results_baselines_covid_t8_full.csv` |
