# Table C.3 — Forecast horizon (paste-ready)

**Status 2026-08-15:** all five jobs **COMPLETED** (7627191–95). Pulled CSVs under `benchmarks/results/table_c3_*_N16.csv`.

**Protocol (locked):** train-H terciles; forced own-H regeneration; Table 5 / Table 4 roots; N=16; K∈{10,100}; methods NeutralBD / ARTreeFormer-adapted / TreeSBM (HIV: Neutral+TreeSBM). See `table_c3_horizon_protocol.md`.

**Coverage@100** below = `coverage_obs_e2` (absolute Hamming ≤2), matching viral main-table convention. CSVs also have e1/e3/e5.

**Tree-KL:** not scored (expensive at K=100) — paper column stays **—** until a dedicated Tree-KL pass.

**Empty buckets (do not fill — locked roots):** COVID short; H1N1 med; H3N2 short; HIV geo med; HIV temporal short.

**Additive layers already in CSV (not in bare paper C.3 header):** `cons_retention`, `aa_acc_given_hit`. **Antigenic** filled from C.3 tree caches (job **7628578** COMPLETED) — see Antigenic column below.

### Locked roots (plain language)

Same Table 5 / Table 4 roots for every horizon. Train-H terciles define short/med/long. Each root is regenerated with its **own** GT mean RTT H, then metrics are averaged inside the tercile that H belongs to. Empty buckets mean none of the locked roots fell there — protocol forbids filling with new roots. Details: `table_c3_horizon_protocol.md`.

---

## Paste-ready (K=100)

| Dataset | Horizon | Method | Cov@100 ↑ | Mut. recall ↑ | Clade recall ↑ | Min dist ↓ | Tree-KL ↓ | Cons ↑ | AA\|hit ↑ | Antigenic ↑ |
|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|
| SARS-CoV-2 Spike | short | — | — | — | — | — | — | — | — | — |
| SARS-CoV-2 Spike | medium | NeutralBD | 0.938 | 0.000 | 0.000 | 0.69 | — | 1.000 | — | 0.010 |
| SARS-CoV-2 Spike | medium | AR | 0.938 | 0.000 | 0.000 | 0.69 | — | 1.000 | — | 0.008 |
| SARS-CoV-2 Spike | medium | TreeSBM | 0.938 | 0.000 | 0.000 | 0.69 | — | 1.000 | — | 0.006 |
| SARS-CoV-2 Spike | long | NeutralBD | 0.667 | 0.003 | 0.000 | 1.90 | — | 1.000 | 1.000† | 0.013 |
| SARS-CoV-2 Spike | long | AR | 0.667 | 0.000 | 0.000 | 1.92 | — | 1.000 | — | 0.005 |
| SARS-CoV-2 Spike | long | TreeSBM | 0.667 | 0.000 | 0.000 | 1.92 | — | 1.000 | — | 0.017 |
| Influenza H1N1 | short | NeutralBD | 1.000 | 0.045 | 0.750 | 0.44 | — | 1.000 | 1.000 | 0.115 |
| Influenza H1N1 | short | AR | 1.000 | 0.067 | 0.250 | 0.46 | — | 1.000 | 1.000 | 0.110 |
| Influenza H1N1 | short | TreeSBM | 1.000 | 0.212 | 0.500 | 0.31 | — | 1.000 | 1.000 | 0.056 |
| Influenza H1N1 | medium | — | — | — | — | — | — | — | — | — |
| Influenza H1N1 | long | NeutralBD | 0.406 | 0.013 | 0.000 | 2.91 | — | 1.000 | 1.000 | 0.107 |
| Influenza H1N1 | long | AR | 0.406 | 0.005 | 0.000 | 2.91 | — | 1.000 | 1.000 | 0.093 |
| Influenza H1N1 | long | TreeSBM | 0.500 | 0.141 | 0.100 | 2.41 | — | 1.000 | 1.000 | 0.097 |
| Influenza H3N2 | short | — | — | — | — | — | — | — | — | — |
| Influenza H3N2 | medium | NeutralBD | 0.828 | 0.000 | 0.000 | 1.27 | — | 1.000 | — | 0.013 |
| Influenza H3N2 | medium | AR | 0.859 | 0.014 | 0.250 | 1.22 | — | 1.000 | 1.000 | 0.026 |
| Influenza H3N2 | medium | TreeSBM | 0.859 | 0.099 | 0.000 | 1.14 | — | 1.000 | 1.000 | 0.011 |
| Influenza H3N2 | long | NeutralBD | 0.000 | 0.000 | 0.000 | 8.25 | — | 1.000 | — | 0.016 |
| Influenza H3N2 | long | AR | 0.000 | 0.000 | 0.000 | 8.25 | — | 1.000 | — | 0.020 |
| Influenza H3N2 | long | TreeSBM | 0.000 | 0.131 | 0.000 | 7.19 | — | 1.000 | 1.000 | 0.005 |
| HIV Env (geo) | short | NeutralBD | 0.313 | 0.001 | 0.067 | 88.6 | — | 1.000 | 0.523 | 0.188 |
| HIV Env (geo) | short | TreeSBM | 0.313 | 0.001 | 0.000 | 88.6 | — | 1.000 | 0.740 | 0.226 |
| HIV Env (geo) | medium | — | — | — | — | — | — | — | — | — |
| HIV Env (geo) | long | NeutralBD | 0.000 | 0.023 | 0.097 | 26.8 | — | 0.982 | 0.576 | 0.142 |
| HIV Env (geo) | long | TreeSBM | 0.000 | 0.050 | 0.097 | 24.7 | — | 0.994 | 0.695 | 0.146 |
| HIV Env (temporal) | short | — | — | — | — | — | — | — | — | — |
| HIV Env (temporal) | medium | NeutralBD | 0.000 | 0.016 | 0.000 | 177.6 | — | 1.000 | 0.694 | 0.212 |
| HIV Env (temporal) | medium | TreeSBM | 0.000 | 0.020 | 0.113 | 177.1 | — | 1.000 | 0.701 | 0.191 |
| HIV Env (temporal) | long | NeutralBD | 0.000 | 0.023 | 0.067 | 186.5 | — | 0.973 | 0.212 | 0.142 |
| HIV Env (temporal) | long | TreeSBM | 0.000 | 0.037 | 0.098 | 184.3 | — | 0.994 | 0.400 | 0.135 |

† sparse / nan-prone when almost no generated mutations hit GT sites.  
Antigenic = mean `lit_hotspot_mut_frac` over gen leaves (PMC / flu lit / HIV V1–V5 masks) from C.3 caches; job **7628578** COMPLETED → `table_c3_*_antigenic.json`.

### Jobs

| Virus | Job | State | Elapsed |
|---|---:|---|---|
| covid | **7627191** | COMPLETED | 05:57:38 |
| h1n1 | **7627192** | COMPLETED | 02:54:29 |
| h3n2 | **7627193** | COMPLETED | 04:26:08 |
| hiv_geo | **7627194** | COMPLETED | 01:30:07 |
| hiv_temporal | **7627195** | COMPLETED | 02:30:36 |
| antigenic score (all) | **7628578** | COMPLETED | 00:03:45 |

### Next (protocol)

1. ~~Pull + aggregate~~ (this file).
2. Tree-KL pass — **deferred** (expensive; optional Day 4–5).
3. Empty-bucket fill — **not allowed** under locked T5/T4 roots.
4. pLM-prior C.3 rows — **not in locked method list**.
5. ~~Antigenic~~ filled from **7628578**.
