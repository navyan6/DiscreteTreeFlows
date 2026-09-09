# Appendix C.1 — Full tree generation (Day 2 paste)

**Updated:** 2026-08-17. HIV Env block replaced by Influenza HA H1N1 (all cells `--` pending **7646375–77**). Do not invent numbers. **—** / `--` = explicit blocker below.  
**Paper columns:** Dataset, Method, Tree-KL↓, Split-KL↓, RF↓, Quartet↓, Branch W1↓, pLM NLL↓  
**Additive (kept on top):** Cons↑, Antigenic↑, AA|hit↑

**Tree-KL honesty:** empirical Tree-KL/Split-KL = NaN in CSVs. Values marked † are **sim_neutral Tree-KL ≈ ln2** (saturated, not discriminative). Prefer RF / Quartet / W1 on empirical; Split-KL on sim_neutral N=16.

Companion paste: [`PASTE_TABLES_2026-08-14.md`](PASTE_TABLES_2026-08-14.md) §3.

---

## Paste-ready

### Protein families × all methods

| Method | Tree-KL | Split-KL | RF | Quartet | W1 | pLM NLL | Cons | Antigenic | AA\|hit |
|---|---|---|---|---|---|---|---|---|---|
| Neutral / AR / TreeSBM | — | — | — | — | — | — | — | — | — |

**Blocker:** no PFAM / protein-family baselines or enrich in this repo.

---

### SARS-CoV-2 Spike — empirical N=16 (`table_empirical_N16.csv` / job **7575017**)

| Method | Tree-KL | Split-KL | RF | Quartet | Branch W1 | pLM NLL | Cons | Antigenic (PMC) | AA\|hit |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Neutral CTMC | 0.693† | 43.0‡ | 0.980 | 0.672 | 2.29e-5 | — | — | — | — |
| Autoregressive (adapted) | 0.693† | 62.9‡ | 0.990 | 0.668 | 5.62e-5 | — | 0.993§ | 0.017§ | 0.094§ |
| TreeSBM | 0.693† | 57.5‡ | 0.982 | 0.669 | 5.45e-5 | **0.457**¶ | **0.839**¶ | **0.005**¶ | **0.375**¶ |

† sim_neutral Tree-KL (saturated). ‡ sim_neutral Split-KL N=16 (`table_sim_neutral_N16.csv`).  
§ Seq from `eval_table5_covid_baselines.json` (job **7539949**) — AR only; Neutral CTMC enrich **not run** (blocker).  
¶ Seq from T8 full enrich `table8_full_mrs0.5.json` / `table7` TreeSBM (geo Brazil mrs=0.5). pLM NLL = 0.4572.

**Neutral CTMC seq (Cons / Ant / aa|hit / NLL):** — · **Blocker:** `eval_table5_covid_baselines.json` has plm_prior + AR only, not Neutral CTMC.

---

### Influenza HA — H3N2 empirical N=16 (`results_baselines.csv` mean over roots)

| Method | Tree-KL | Split-KL | RF | Quartet | Branch W1 | pLM NLL | Cons | Antigenic (flu lit) | AA\|hit |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Neutral CTMC | 0.693† | 43.0‡ | 0.980 | 0.714 | 7.77e-4 | **0.466** | — | — | — |
| Autoregressive (adapted) | 0.693† | 62.9‡ | 0.986 | 0.717 | 9.21e-4 | **0.462** | — | — | — |
| TreeSBM | 0.693† | 57.4‡ | 0.983 | 0.721 | 9.25e-4 | **0.463** | **0.918**¶ | **0.006**¶ | **0.439**¶ |

†‡ sim_neutral N=16 means from same CSV (`split_kl` ≈ 42.97 / 62.89 / 57.36).  
¶ `eval_enrichment_h3n2_v3_lit_hotspot_mrs0.5_fullmetrics.json` (`lit_hotspot_mut_frac` = flu primary sites).  
**pLM NLL:** Neutral/AR from `eval_table5_h3n2_baselines_plmnll.json`; TreeSBM from **7635752** `eval_enrichment_h3n2_v3_lit_hotspot_mrs0.5_plmnll.json` (**0.4627**).  
**Neutral / AR Cons/Ant/aa|hit:** — · **Blocker:** no H3 Neutral/AR enrich JSON with those fields (baselines_plmnll has cons/aa|hit but not lit ant).

---

### Influenza HA H1N1 — geo empirical N=16 (replaces HIV Env in C.1 `.tex`)

**Best TreeSBM:** `checkpoints/h1n1_v2_lit_hotspot/best.pt` (ep74) on **geo** `data/h1n1` — same ckpt/split as Table 5 / C.3. **Not** `h1n1_leafholdout` / `h1n1_temporal` (those exist but were unused for paper T5).

| Method | Tree-KL | Split-KL | RF | Quartet | Branch W1 | pLM NLL | Cons | Antigenic (H1 lit) | AA\|hit |
|---|---|---|---|---|---|---|---|---|---|
| Neutral CTMC | -- | -- | -- | -- | -- | -- | -- | -- | -- |
| Autoregressive (adapted) | -- | -- | -- | -- | -- | -- | -- | -- | -- |
| TreeSBM (`h1n1_v2_lit_hotspot`) | -- | -- | -- | -- | -- | -- | -- | -- | -- |

**Do not invent.** No prior H1 topology suite (`results_baselines_h1n1.csv` absent). `eval_table5_h1n1_baselines.json` (**7539950**) is pLM+AR only and has **no** `plm_nll`. TreeSBM mrs=0.5 enrich has **no** `plm_nll`.

**Queued 2026-08-17 (Betty login02):** [`table_c1_h1n1_job_ids.json`](table_c1_h1n1_job_ids.json)

| Job | Name | Fills | Artifact |
|---:|---|---|---|
| **7646375** | `c1h1_bl` | Tree-KL / Split-KL / RF / Quartet / Branch W1 (all 3 methods) | `results_baselines_h1n1.csv` + `tables/h1n1/table_empirical_N16.csv` |
| **7646376** | `c1nll_h1` | Neutral CTMC + AR pLM NLL | `checkpoints/eval_table5_h1n1_baselines_plmnll.json` |
| **7646377** | `c1nll_h1ts` | TreeSBM pLM NLL | `checkpoints/eval_enrichment_h1n1_v2_lit_hotspot_mrs0.5_plmnll.json` |

Analogous to C.1 NLL **7635749–52** (covid/h3/hiv) + topo **7635753** / **7626610**. Pull when COMPLETED; leave `--` until then.

HIV Env C.1 numbers are **not** in the LaTeX table anymore; they remain in `hiv_{geo,temporal}/` + [`hiv_c1_job_ids.json`](hiv_c1_job_ids.json).

---

### Antibody lineages

| Method | Tree-KL | Split-KL | RF | Quartet | W1 | pLM NLL | Cons | Antigenic | AA\|hit |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| Neutral CTMC | **0.693** | **113.1** | **0.970** | **0.639** | **0.00670** | **0.513** | — | — | — |
| Autoregressive (adapted) | **0.693** | **115.4** | **0.984** | **0.648** | **0.00859** | **0.466** | — | — | — |
| TreeSBM (OAS `ab_oas_1m_v1`) | **0.693** | **118.1** | **0.979** | **0.658** | **0.0109** | **0.451** | — | — | — |

**Source:** free-topo job **7635753** → `results_baselines_ab_oas.csv` / `tables/ab_oas/table_empirical_N16.csv` (n_roots=12). NLL from **7635754/55**. Track C / Rod.82 remains forced-topo and is **not** the C.1 topo source.

---

## Artifact map

| Dataset | Topo | Seq enrich | Notes |
|---|---|---|---|
| COVID | `results_baselines_covid.csv` **7575017**; agg `table_empirical_N16.csv` | T8/T7 TreeSBM; `eval_table5_covid_baselines.json` AR/pLM | Neutral seq missing |
| H3N2 | `results_baselines.csv` (mean N=16 emp) | `eval_enrichment_h3n2_*_mrs0.5_fullmetrics.json` | Neutral/AR seq missing |
| H1N1 geo | **7646375** `results_baselines_h1n1.csv` (PENDING) | NLL **7646376/77** PENDING · ckpt `h1n1_v2_lit_hotspot` | all C.1 cells `--` until COMPLETED |
| HIV geo/temporal | `hiv_*/table_empirical_N16.csv` | `eval_table5_hiv_*_baselines.json` + TreeSBM mrs0.5 | COMPLETED but **removed from C.1 `.tex`** |
| Ab | free-topo **7635753** | NLL **7635754/55** | Track C forced topo unused |
| Protein families | — | — | no data |
