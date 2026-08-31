# EVEREST clade forecasting eval

Evaluates TreeSBM (generative mode) against the [priority-viruses forecasting_dataset](https://github.com/debbiemarkslab/priority-viruses/tree/main/data/forecasting_dataset). **No** `prepare_everest_clades.py` — CSVs are read inline at eval time.

## Script

`scripts/eval_everest_forecasting.py`

Clone data once:

```bash
git clone --depth 1 https://github.com/debbiemarkslab/priority-viruses external/priority-viruses
```

## Metrics (per clade CSV × checkpoint)

| Metric | Definition |
|--------|------------|
| `emerging_mut_recall@K` | Fraction of clade mutations hit by ≥1 generated leaf (correct AA, root≠alt) |
| `weighted_recall` | Same, weighted by `count` (flu) or `test_count` (COVID/HIV) |
| `enrichment_vs_random` | `weighted_recall` / mean permuted-position null (keeps mutant AA labels) |
| `mut_recovery` | Best-match GT leaf vs root (standard TreeSBM leaf metric) |
| `coverage_obs_e2` | Fraction of GT leaves within Hamming ≤2 of some gen leaf |

## Coordinate mapping

| Virus | CSV mutation | Column index |
|-------|--------------|--------------|
| H3N2 / H1N1 HA | `I10T` (H3/H1 mature 1-based) | `col = pos + signal_len - 1` (H3 signal=16, H1 default=17) |
| SARS-CoV-2 Spike | `A1015D` (Wuhan 1-based) | `col = pos - 1` |
| HIV Env | `M1K` | `col = pos - 1` |

Clade window: filename suffix `_YYYY-MM-DD_YYYY-MM-DD`; roots with collection date after window end are skipped. Mutations with `first_seen` after window end are filtered out.

## Checkpoint ↔ data ↔ EVEREST folder

| Checkpoint | Test data | EVEREST folder |
|------------|-----------|----------------|
| `covid_v6_epidemic_mutrec` | `data/covid_epidemic/test` | `SARS-CoV-2/` |
| `h3n2_v4_epidemic_mutrec` | `data/h3n2_epidemic/test` | `H3N2/` |
| `h1n1_v3_epidemic_mutrec` | `data/h1n1_epidemic/test` | `H1N1/` |
| Wave-1 ablation ckpts | geo/temporal test dirs | same folders |

## Betty batch

```bash
bash scripts/slurm_eval_everest.sh
# or sbatch wrapper if added
```

Outputs: `benchmarks/results/everest/<name>.json`

## Pan-flu forecast generality

After pan-flu holdout trains (`betty_submit_panflu_forecast.sh`), run EVEREST on H3N2/H1N1 clades with 2020 `first_seen` and compare B0 (H3N2-only train) vs G3 (pan-flu train) on the same H3N2 2020 test.

See [`SPLITS.md`](SPLITS.md) § Pan-flu forecast generality.
