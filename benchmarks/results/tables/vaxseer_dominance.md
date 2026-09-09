# VaxSeer dominance validation (flu HA only)

Status: **repo cloned + official results dump extracted**; live **lm weight download blocked** on interactive Dropbox token.

## Setup (Betty)

| Piece | Path |
|---|---|
| Clone | `$LABHOME/third_party/vaxseer` ([wxsh1213/vaxseer](https://github.com/wxsh1213/vaxseer.git)) |
| Paper results dump | `$LABHOME/third_party/vaxseer_results/runs/pipeline/{Y}-02/{a_h3n2\|a_h1n1}/…/dominance_prediction/lm/…/test_results.csv` |
| Adapter | `benchmarks/methods/vaxseer_dominance.py` |
| Eval | `scripts/eval_vaxseer_dominance.py` |
| Setup helper | `scripts/setup_vaxseer.sh` |

### lm weights (required for scoring *generated* TreeSBM leaves)

```bash
cd $LABHOME/third_party/vaxseer
# needs Dropbox access token (interactive prompt in download_models_from_dropbox.py)
python download_models_from_dropbox.py --task lm --year 2018 --subtype a_h3n2 --output_dir runs
```

Years available upstream: **2012–2021** only (not 2024). Align TreeSBM temporal H3 cutoff to a VaxSeer year (e.g. 2018/2019) or regenerate from matching pre-Feb roots.

Precomputed `test_results.csv` rows are keyed by **GISAID EPI ids**, not AA strings — useful for observed panels, not for novel generated leaves. Live LM ckpt required for generated \(p_Y(g)\).

## Protocol

For each subtype ∈ {H3N2, H1N1} and season \(Y\) (train cutoff ≈ before 1 Feb of \(Y\)):

| Set | Definition | Dominance |
|---|---|---|
| Observed \(O_Y\) | Circulating HA in season \(Y\) | Empirical \(f_Y\) and VaxSeer \(p_Y\) |
| Generated \(G\) | TreeSBM / Neutral / pLM / AR leaves from pre-season roots | VaxSeer \(p_Y\) only |

Primary reports: top-decile recovery @ε∈{2,5}; \(p_Y\) mean/median/max; calibration when freqs available.

COVID/HIV: skip (HA-specific).

## Results

| subtype | year | lm weights | TreeSBM top10%@ε2 | Neutral | pLM | AR | notes |
|---|---:|---|---:|---:|---:|---:|---|
| a_h3n2 | 2018 | pending Dropbox token | — | — | — | — | results dump has EPI-id preds |
| a_h1n1 | 2018 | pending Dropbox token | — | — | — | — | same |

Do not invent numbers. After token download: `eval_vaxseer_dominance.py` → fill this table + PASTE/APPENDIX.
