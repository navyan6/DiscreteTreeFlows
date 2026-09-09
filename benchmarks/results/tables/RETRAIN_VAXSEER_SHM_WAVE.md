# Retrain + VaxSeer + Ab SHM wave — numbers tracker

Update when jobs complete. **Do not invent numbers. Do not overwrite paper ckpts.**

Submitted Betty **2026-08-20** from login01 (`betty_submit_retrain_wave.sh`).

## Locked paper ckpts (never overwrite)

`covid_v5_mutrec` · `h1n1_v2_lit_hotspot` · `h3n2_v3_lit_hotspot` · `hiv_geo_v1` · `ab_oas_1m_v1` · pathogen `best.pt`

## Viral eval-only sweep

Job IDs: [`viral_eval_sweep_job_ids.json`](viral_eval_sweep_job_ids.json) (**7720891–7720910**)

Paper-ckpt mrs=0.5 enrich (from `WAVE_2026-08-22_RESULTS.md` §1): COVID **0.069 / 0.341**; H1 **0.154 / 0.618**; H3 **0.100 / 0.439**; HIV **0.154 / 0.310** (mut / aa\|hit). Clade_recall is coverage-CSV only, not enrich.

Fill from `checkpoints/viral_eval_sweep/*_enrich.json` and `benchmarks/results/viral_eval_sweep_*.csv`.

## Viral retrain mutlin (β∈{0,0.25})

Job IDs: [`viral_retrain_mutlin_job_ids.json`](viral_retrain_mutlin_job_ids.json) (**7720911–7720918**)

Enrich mrs=0.5 jobs **7794867–7794874** COMPLETED. Compact: [`viral_mutlin_enrich_summary.json`](viral_mutlin_enrich_summary.json). No clade_recall in enrich JSONs.

| ckpt | β | job | val_loss | mut_recovery | aa\|hit |
|---|---:|---:|---:|---:|---:|
| covid_v6_mutlin_b0 | 0 | 7720911 | 0.651 | 0.0698 | 0.277 |
| covid_v6_mutlin_b025 | 0.25 | 7720915 | 0.810 | 0.0554 | 0.306 |
| h1n1_v3_mutlin_b0 | 0 | 7720912 | 3.551 | 0.124 | 0.558 |
| h1n1_v3_mutlin_b025 | 0.25 | 7720916 | 3.777 | 0.139 | 0.628 |
| h3n2_v4_mutlin_b0 | 0 | 7720913 | 5.151 | 0.118 | 0.416 |
| h3n2_v4_mutlin_b025 | 0.25 | 7720917 | 6.087 | 0.182 | 0.685 |
| hiv_geo_v2_mutlin_b0 | 0 | 7720914 | 76.67 | 0.166 | 0.359 |
| hiv_geo_v2_mutlin_b025 | 0.25 | 7720918 | 79.12 | 0.149 | 0.348 |

## Flu variant path

- Tree viz **7720921** → `results/h3n2_tree_viz/`
- Path figure afterok **7720924** → `results/h3n2_variant_path/`
- Script: `scripts/build_flu_variant_path_figure.py`

## VaxSeer

See [`vaxseer_dominance.md`](vaxseer_dominance.md). Clone + results dump done; **lm Dropbox token** still needed for scoring generated leaves.

## Ab SHM

- Entropy audit **7720923** COMPLETED — [`ab_cdr_fwr_entropy_audit.md`](ab_cdr_fwr_entropy_audit.md)
  - OAS H_CDR/H_FWR = **1.627**; Rod.82 = **2.024** (CDR-biased; supports Recipe A)
- v2 train **7720919** → `checkpoints/ab_oas_1m_v2_shm`
- Rod.82 **7720920** → `treesbm_ab_oas_v2`
- IDs: [`ab_oas_v2_shm_job_ids.json`](ab_oas_v2_shm_job_ids.json)
- Recipe B **COMPLETED** 7799979 / 7799980 / 7799981. Ckpt `ab_oas_1m_v3_thrifty` (val 71.33). Fair T6: [`table6_ab_oas_v2v3.md`](table6_ab_oas_v2v3.md).
- Recipe C **COMPLETED** 7816338 / 7816339. Ckpt `ab_oas_1m_v4_cosinehead` (val 134.73). Fair T6: [`table6_ab_oas_v2v3v4.md`](table6_ab_oas_v2v3v4.md). CDR **0.158**.

| model | CDR recall | SHM load err | term-div err |
|---|---:|---:|---:|
| treesbm_ab_oas (v1, N=20) | 0.150 | 0.077 | 9.89 |
| treesbm_ab_oas_v2 (Recipe A, N=20) | 0.162 | 0.078 | 9.98 |
| treesbm_ab_oas_v3 (Thrifty Q0, N=20) | **0.175** | 0.083 | 11.13 |
| treesbm_ab_oas_v4 (cosine head, N=20) | 0.158 | 0.054 | 4.89 |
| Neutral JC69 | 0.688 | 0.054 | 17.44 |
| CoSiNE | 0.658 | 0.037 | 3.02 |
| Thrifty | 0.746 | 0.037 | 12.97 |

## PASTE / APPENDIX

When numbers land: copy into `PASTE_TABLES_*.md` and `APPENDIX_COMPLETION_PLAN.md` (Standing SOP: update **numbers**). Tracker for this wave: this file.
