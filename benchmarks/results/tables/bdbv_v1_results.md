# BDBV L TreeSBM v1 results (scaffold)

Fill after train + eval complete. Job IDs: `bdbv_job_ids.json`.

## Primary metrics (temporal test)

| Checkpoint | Cov@K (ε=5) | mut_recovery | site_recall | aa\|hit | cons | K1738N aa\|hit | Q1770R aa\|hit |
|------------|------------:|-------------:|------------:|--------:|-----:|---------------:|---------------:|
| bdbv_v1_mutrec | — | — | — | — | — | — | — |
| bdbv_v1_lit_mutrec | — | — | — | — | — | — | — |
| bdbv_pan_v1_mutrec | — | — | — | — | — | — | — |
| bdbv_pan_v1_lit_mutrec | — | — | — | — | — | — | — |

## Lit stratification (v1 vs v1_lit)

| Checkpoint | lit mut_recovery | bg mut_recovery | lit site_recall | bg site_recall |
|------------|-----------------:|----------------:|----------------:|---------------:|
| bdbv_v1_mutrec | — | — | — | — |
| bdbv_v1_lit_mutrec | — | — | — | — |

## Decision table

| Question | Result |
|----------|--------|
| Temporal vs geo | pending SPLIT_AUDIT |
| Pan vs BDBV-only | pending |
| Lit mask on vs off | pending (cons ≥ 0.98 gate) |

## Artifacts

- Coverage: `benchmarks/results/coverage_curves_bdbv_N16_eabs.csv`
- Enrichment: `checkpoints/eval_enrichment_bdbv_v1_mutrec_mrs*.json`
- Variants: `checkpoints/eval_bdbv_2026_variants_bdbv_v1_mutrec.json`
- Fixed topo: `checkpoints/viral_eval_sweep/bdbv_fixed_topo_mrs1.0.json`
