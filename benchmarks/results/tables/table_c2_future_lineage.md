# Appendix C.2 — Future-lineage coverage (`tab:app_future_lineage_full`)

**Updated:** 2026-08-25 (COVID TreeSBM Cov@500 filled from `coverage_curves_covid_N16_table5_eabs_K500.csv` `coverage_obs_e2` @ K=500 = **0.775**; prior job **7635911**/CSV now on disk).  
**Paper columns:** Dataset, Method, Cov@100↑, Cov@500↑, Mut. recall↑, Clade recall↑, Min dist↓  
**LaTeX:** [`table_c2_future_lineage.tex`](table_c2_future_lineage.tex) — structure unchanged; do not invent numbers.  
**Naming:** User “c.3” paste = this table (Appendix **C.2**). Horizon = **C.3** (`table_c3_horizon.md`).

**Cov@100** = `coverage_obs_e2` at K=100 (absolute Hamming ≤ 2) unless noted.  
**Cov@500** = same metric at K=500 from dedicated K500 curves (not fractional ε).

Mut. recall prefers **enrich** JSONs when available; Clade recall + Min dist from **coverage** CSVs (K=100). Ab Mut. recall = **CDR recall** (Track C). Ab Cov/Clade/Min from `table_c2_ab_metrics.json` (**7639280**).

---

## Paste-ready (paper 5 metric cols)

| Dataset | Method | Cov@100 | Cov@500 | Mut | Clade | Min dist |
|---|---|---:|---:|---:|---:|---:|
| SARS-CoV-2 Spike | pLM prior | 0.775 | **0.775** | 0.048 | 0.000 | 1.425 |
| SARS-CoV-2 Spike | TreeSBM | 0.788 | **0.775** | 0.080 | 0.000 | 1.388 |
| Influenza HA H3N2 | pLM prior | 0.675 | **0.750** | 0.004 | 0.000 | 2.650 |
| Influenza HA H3N2 | TreeSBM | 0.725 | **0.788** | 0.100 | 0.100 | 2.238 |
| HIV Env | pLM prior | **0.188** | **0.188** | **0.015** | **0.064** | **63.9** |
| HIV Env | TreeSBM | 0.188 | **0.200** | 0.026 | 0.061 | 63.1 |
| Antibody lineages | pLM prior | 0.018 | **0.018** | 0.297 | **0.039** | **15.9** |
| Antibody lineages | TreeSBM | **0.027** | **0.027** | 0.150 | **0.018** | **14.3** |

---

## Source picks (best TreeSBM / pLM)

| Dataset | TreeSBM pick | pLM | Notes |
|---|---|---|---|
| SARS-CoV-2 | `covid_v5_mutrec` (ep14) | Table 5 enrich + cov **7575015** / K500 CSV | Cov@100 from `coverage_curves_covid_N16_table5_eabs.csv`; Mut from enrich. Cov@500 TreeSBM = **0.775** (`coverage_curves_covid_N16_table5_eabs_K500.csv` `coverage_obs_e2` @ K=500; same plateau as pLM — K does not lift mid-ε COVID) |
| H3N2 | `h3n2_v3_lit_hotspot` | cov `coverage_curves_h3n2_N16_eabs.csv` | Cov@100 from original eabs K=100. Mut TreeSBM = enrich mrs0.5 **0.100**. Cov@500 from rescore of K500 cache |
| HIV Env | **geo** `hiv_geo_v1` | **7635912–14** COMPLETED | Cov@100/500 e2=**0.188** (`…_eabs_plm.csv` / `…_K500_plm.csv`). Clade/Min @K=100. Mut=`mut_recovery` @K=100 cov (**0.015**; enrich JSON lacks mut_recovery) |
| Antibody | **`treesbm_ab_oas`** | Rod.82; **7639280** | Cov/Clade/Min from `table_c2_ab_metrics.json`. Cov@100 OAS **0.027** (100 rollouts; was 0.024 @ N=20 Track C). Mut = CDR (Track C unchanged) |

### Additive (not in bare paper tex; for PASTE)

| Dataset | Method | Cons | Antigenic | AA\|hit |
|---|---|---:|---:|---:|
| COVID | pLM / TreeSBM | 0.985 / 0.839 | 0.000 / 0.005 | 0.017 / 0.375 |
| H3N2 | pLM / TreeSBM | 1.000† / **0.918** | — / **0.006** | — / **0.439** |
| HIV geo | pLM / TreeSBM | 0.770 / 0.893 | 0.143 / 0.134 | 0.077 / 0.337 |
| Ab | pLM / OAS TreeSBM | — | — | — |

† H3 pLM cons from coverage CSV (no enrich).

---

## Remaining `--` + job IDs

| Cell | Status | Job / action |
|---|---|---|
| COVID TreeSBM **Cov@500** | — → RUNNING (~5.2h) | **7635911** on dgx028; NODE_0000179 = Table5 group_010 (present in test); treesbm cache still 4/5 (write-after-K_max). Not FAILED/stuck — leave running; Cov CSV still empty |
| HIV pLM all cols | **FILLED** | **7635912–14** COMPLETED |
| Ab Cov@500 / Clade / Min (pLM + OAS) | **FILLED** | rolls **7635915/16** + eval **7639280** (comma-export bug on **7635917** fixed) |

### Cov@500 wave (**7628571–75**) + gap resume

| Virus | Job | Status | Artifact |
|---|---:|---|---|
| COVID | **7628571** → **7635911** | cancelled incomplete; **RUNNING** resume | `coverage_curves_covid_N16_table5_eabs_K500.csv` |
| H1N1 | **7628572** | DONE (not in paper C.2 rows) | `coverage_curves_h1n1_N16_table5_eabs_K500.csv` |
| H3N2 | **7628573** | DONE gen; rescored e-abs | Cov@500 pLM=0.750 TreeSBM=0.788 |
| HIV geo TreeSBM | **7628574** | DONE | Cov@500=0.200 |
| HIV geo pLM | **7635912/13/14** | **COMPLETED** | Cov@100=Cov@500=0.188; clade=0.064; min=63.9; mut=0.015 |
| HIV temporal | **7628575** | DONE (geo preferred) | — |
| Ab | **7635915–17** / **7639280** | **COMPLETED** | `table_c2_ab_metrics.json` |

Job ID JSON: [`table_c2_cov500_job_ids.json`](table_c2_cov500_job_ids.json).

---

## C.3 note

Horizon table is Appendix **C.3** (`table_c3_horizon.md`) — already filled for non-empty buckets. This file is **C.2** future-lineage only.
