# Table 5 — Viral forecasting (2026-08-12 e-abs pull)

**Coverage = `coverage_obs_e2` (absolute Hamming ≤ 2).** Fractional ε=2% saturates at 1.0 — do not use.

**COVID pLM EVEscape Δ filled 2026-08-25:** `checkpoints/eval_table5_covid_plm_evescape.json` → `evescape_score_delta` = **−0.424** (job **7644883** path; L=1280 fix).

## Paste-ready

| Virus | Method | Cov@10 | Cov@100 | Key mut | Antigenic | AA\|hit | Future rank | Min dist | EVEScape Δ | Cons |
|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
| SARS-CoV-2 Spike | pLM mutation prior | 0.775 | 0.775 | 0.048 | 0.000 | 0.017 | — | 1.425 | **−0.424** | 0.985 |
| SARS-CoV-2 Spike | Autoregressive tree-edit | 0.775 | 0.775 | 0.006 | 0.017 | 0.094 | — | 1.425 | −0.371 | 0.993 |
| SARS-CoV-2 Spike | TreeSBM | 0.775 | 0.788 | 0.080 | 0.005 | 0.375 | — | 1.388 | −0.048 | 0.839 |
| Influenza H1N1 HA | pLM mutation prior | 0.762 | 0.762 | 0.071 | 0.115 | 0.136 | — | 1.462 | −0.354 | 0.945 |
| Influenza H1N1 HA | Autoregressive tree-edit | 0.762 | 0.762 | 0.006 | 0.000 | 0.008 | — | 1.462 | −0.298 | 0.998 |
| Influenza H1N1 HA | TreeSBM | 0.800 | 0.762 | 0.241 | 0.163 | 0.618 | — | 1.388 | −0.138 | 0.952 |

## Protocol

| | COVID | H1N1 |
|---|---|---|
| Coverage jobs | **7575015** COMPLETED | **7575016** COMPLETED |
| Ckpt | `covid_v5_mutrec` (ep14) | `h1n1_v2_lit_hotspot` (ep74) |
| Split | Brazil geo-test | location-unit geo-test |
| Coverage roots (N=16, 5) | groups **004,005,007,009,010** (Brazil, 2021-05-10…2021-07) | **001** Canada:Ontario; **004** Iran; **005** Iran; **006** Canada-pooled; **007** Canada:Toronto |
| Enrichment | 7539949 + 7539980 (groups 1–20); COVID pLM EVEscape **7644883** → `eval_table5_covid_plm_evescape.json` | 7539950 + prior h1n1_v2 mrs=0.5 |
| L | 1280 | 566 |
| EVEscape Δ | `model_evescape − gt_evescape` (all scored muts) | same |

CSV: `benchmarks/results/tables/table5_viral_forecasting.csv`
