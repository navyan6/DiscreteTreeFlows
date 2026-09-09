# Table 4: Future lineage forecasting (ε = 5 AA Hamming)

**Setting:** H3N2 N=16, 5 held-out roots, job **7459412**, artifact `benchmarks/results/coverage_curves_h3n2_N16_eabs.csv`.

**Coverage** = `coverage_obs_e5` (fraction of observed leaves within Hamming ≤ 5 of some generated leaf).

**Mutation recall** = `mut_recovery` from the same coverage run at **K=100** (`aa_acc_given_hit` not emitted by coverage_curves).

**Clade recall** / **min dist** = `clade_recall` / `mean_min_edit` at **K=100**. Paper column says “Median min dist”; code reports **mean** min Hamming over observed leaves.

| Method | Coverage@10 ↑ | Coverage@100 ↑ | Mutation recall ↑ | Clade recall ↑ | Median min dist. ↓ |
|---|---:|---:|---:|---:|---:|
| Neutral CTMC + birth-death | 0.800 | 0.800 | 0.020 | 0.000 | 2.65 |
| pLM mutation prior only | 0.800 | 0.800 | 0.004 | 0.000 | 2.65 |
| Autoregressive tree-edit model | 0.800 | 0.812 | 0.079 | 0.300 | 2.40 |
| TreeSBM without bridge control | — | — | — | — | — |
| TreeSBM | 0.800 | 0.812 | 0.156 | 0.100 | 2.24 |

## Method mapping

| Paper name | Code `method` |
|---|---|
| Neutral CTMC + birth-death | `neutral_bd` |
| pLM mutation prior only | `plm_prior` |
| Autoregressive tree-edit model | `artreeformer_adapted` |
| TreeSBM without bridge control | — (no run) |
| TreeSBM | `treesbm` (`h3n2_v3_lit_hotspot`) |

## Gaps

- **TreeSBM without bridge control:** no e-abs coverage artifact on Betty/local (`reference_rollout` exists as a script but was not run in this protocol).
- Paper “Median min dist” vs code `mean_min_edit` (mean, not median).
- `aa_acc_given_hit` unavailable in coverage CSV; enrichment JSONs use a different protocol (not N16 forecasting).

## Plot series (no plots generated)

- `table4_coverage_vs_eps.csv` — coverage vs Hamming ε at K∈{10,100}
- `table4_coverage_vs_K_eps5.csv` — coverage_obs_e5 / mut_recovery / clade_recall / mean_min_edit vs K
- `table4_coverage_by_eps_wide.csv` — wide form for quick paste

