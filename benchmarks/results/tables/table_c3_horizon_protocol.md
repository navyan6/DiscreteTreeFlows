# Table C.3 — genetic horizon stratification (protocol note)

**Design (locked):** Horizon = train-set genetic H terciles (`h_buckets`). Regenerate with forced H = each example’s own mean RTT (true conditioning). Stratify metrics by which tercile that H falls in. Same Table 5 roots (COVID / H1 / HIV) + Table 4 H3N2 roots. N=16, K∈{10,100}.

## Locked roots (plain language)

We do **not** pick new roots per horizon. We reuse the **same five Table 5 / Table 4 evaluation roots** already used for main-text coverage:

1. Compute train-set genetic depth H (mean root-to-tip) tercile cuts → buckets **short / med / long**.
2. For each locked root, read its **own** GT H and assign it to whichever bucket that H falls in.
3. When regenerating, pass that **same own H** into `generate(root_seq, N, H)` (true conditioning), then report metrics **stratified by bucket**.

So “locked” means: no expanding empty buckets with extra roots; empty cells (e.g. COVID short) stay empty because **none of the five fixed roots** land in that tercile. Horizon is a property of the root’s depth, not a free sampling knob.

## Pathogens / ckpts / L

| Virus | Roots | Ckpt | max_seq_len |
|---|---|---|---|
| COVID Brazil geo | groups 4,5,7,9,10 | `covid_v5_mutrec` | 1280 |
| H1N1 geo | groups 1,4,5,6,7 | `h1n1_v2_lit_hotspot` | 566 |
| H3N2 temporal (primary flu) | group 2 ×5 (Table 4) | `h3n2_v3_lit_hotspot` | 566 |
| HIV geo (prefer) | diverse 5 (T5 filter) | `hiv_geo_v1` | 900 |
| HIV temporal | diverse 5 | `hiv_temporal_v1` | 900 |

## Methods / metrics

- **Methods:** NeutralBD, ARTreeFormer-adapted (COVID/flu; HIV: Neutral+TreeSBM), TreeSBM
- **Metrics:** Coverage@K `e1/e2/e3/e5`, `mut_recovery`, `clade_recall`, `mean_min_edit`
- **Tree-KL:** not run (expensive at K=100); add later if needed
- **Outputs:** `table_c3_{virus}_N16.csv` (by bucket), `*_per_root.csv`, `table_c3_{virus}_h_cuts.json`

## Submit

```bash
VIRUS=covid sbatch --qos=mig-max scripts/slurm_table_c3_horizon.sh
# … h1n1, h3n2, hiv_geo, hiv_temporal
```

Job IDs: `benchmarks/results/tables/table_c3_horizon_job_ids.json`
