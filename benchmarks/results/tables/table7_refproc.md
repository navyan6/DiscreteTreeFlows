# Table 7 — Reference process / R0 backends (COVID Brazil geo-test)

**Updated 2026-08-16** (ESM-C / ProGen2 Cov@100 **7628260/61/62** COMPLETED — filled Cov e1/e2/e5 + min_edit). Paste companion: [`PASTE_TABLES_2026-08-14.md`](PASTE_TABLES_2026-08-14.md).
Data `data/covid/test`; ckpt `checkpoints/covid_v5_mutrec/best.pt`; mrs=0.5; max_trees=20; L=1280.
Rows 1–6: `--ablate-bridge` (pure R0). Row 7: full TreeSBM + `--fitness-beta 1.0`.

**Coverage roots (all T7/T8 Cov jobs):** Brazil groups **4, 5, 7, 9, 10** — NODE_0000162 / 282 / 271 / 0000000 / 0179. Group 9 mean_edit=5.62 drives Cov@e2=0 on that root; groups 5/7/10 are near-identical (edit≈0.06–0.31) so **Cov@100 e2 ≈ 0.775 for almost every method** — weak discriminator. Prefer mut-recall / aa|hit / min_edit. ESM-C / ProGen2 Cov CSVs: `coverage_curves_covid_N16_table7_{esmc_nofit,esmc_fit,progen2_fit}_K100.csv`.

**Tree-KL / Split-KL / W1:** empirical track is NaN. Values below for TreeSBM are **7575017 COMPLETED** `sim_neutral` N=16 (n=20 roots). Not computed per R0 backend (T7 coverage jobs do not emit KL). Tree-KL ≈ ln2 for *all* methods in 7575017 — saturated.

**pLM NLL:** JTT **7626622**=0.4705; ESM2-nofit **7626623**=0.4415; ESM2-fit **7626624**=0.4405; TreeSBM **7626625**=0.4464; ESM-C-nofit **7627569**=0.4605; ESM-C-fit **7627570**=0.4267; ProGen2 **7627918**=0.4546 (**n=20**; supersedes partial **7627571** n=3).

| R0 / method | Enrich job | Cov job | State | Mut recall | site_recall | aa\|hit | Cons | Antigenic (PMC) | EVE Δ | Cov@100 e1 | **Cov@100 e2** | Cov@100 e5 | min_edit@100 | Tree-KL | Split-KL | Branch W1 | pLM NLL |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| JTT/WAG/LG (no fitness) | 7626622 | 7585984 | COMPLETED | 0.0223 | 0.2330 | 0.0825 | 0.8213 | 0.0066 | -0.0702 | 0.7750 | **0.7750** | 0.9625 | 1.4125 | — | — | — | **0.4705** |
| ESM-2-650M (no fitness) | 7626623 | 7585985 | COMPLETED | 0.0180 | 0.1640 | 0.1083 | 0.8931 | 0.0059 | -0.0808 | 0.7750 | **0.7750** | 0.9750 | 1.4125 | — | — | — | **0.4415** |
| ESM-2-650M (fitness) | 7626624 | 7585986 | COMPLETED | 0.0011 | 0.0196 | 0.0983 | 0.9884 | 0.0018 | -0.4553 | 0.7750 | **0.7750** | 0.9625 | 1.4250 | — | — | — | **0.4405** |
| ESM-C (no fitness) | **7627569** | **7628260** | COMPLETED | 0.0347 | 0.3631 | 0.0796 | 0.6923 | 0.0055 | -0.0370 | 0.7750 | **0.7750** | 0.9750 | 1.2250 | — | — | — | **0.4605** |
| ESM-C (fitness) | **7627570** | **7628261** | COMPLETED | 0.0356 | 0.3915 | 0.0853 | 0.7134 | 0.0053 | -0.0332 | 0.7750 | **0.7750** | 0.9625 | 1.4125 | — | — | — | **0.4267** |
| ProGen2 (fitness) | **7627918** | **7628262** | COMPLETED n=20 | 0.1126 | 0.1365 | 0.3978 | 0.9855 | 0.0048 | -0.0325 | 0.7750 | **0.7875** | 0.9750 | 1.1625 | — | — | — | **0.4546** |
| TreeSBM default (fitness β=1) | 7626625 | 7585987 (repeat of 7575916) | COMPLETED | 0.0378 | 0.1047 | 0.2896 | 0.9412 | 0.0061 | -0.0936 | 0.7750 | **0.7750** | 0.9750 | 1.2375 | 0.6931 | 57.5201 | 0.000052 | **0.4464** |

‡ **ProGen2 prior (7627571):** COMPLETED but only **3/20 trees** — **superseded** by **7627918** (n=20 after chunked long-seq fix + smoke **7627917**). Do **not** paste 7627571 as primary.

**R0 backend fix (2026-08-14 / length 2026-08-15)**

- **ESM-C:** `src/r0_backends.py` now loads HF `Synthyra/ESMplusplus_small` via `AutoModelForMaskedLM` (`attn_implementation=sdpa`/`eager`). No fair-esm; slurm no longer sets `PYTHONPATH=python_esmc`. Official `biohub/ESMC-300M` still needs native transformers `esmc` (not in 5.13–5.15).
- **ProGen2:** official [enijkamp/progen2](https://github.com/enijkamp/progen2) at `/vast/projects/pranam/lab/nnori/progen2` + GCS `progen2-small`. Likelihood-style next-token AA log-probs. transformers≥5 shims for `get_head_mask` / `model_parallel_utils`; weights via `load_state_dict` (not hugohrban AutoModel).
- **Long Spike (2026-08-15):** `PROGEN2_MAX_CTX=1024`; sequences longer than 1023 AA use **sliding left-context chunk aggregation** (not silent truncate). Prefer `progen2-small`. Smoke job **7627917** (L=1273 → shape `(1,1280,20)`); enrich **7627918**.
- Smoke (Betty MIG): short-seq both backends → `(1, 32, 20)`; long-seq ProGen2 → **7627917**.
- Historical failures: **7585981/82** libcudart; **7585983** `get_head_mask` / 0 trees; **7627571** partial n=3.

**7575017** COVID all-method baselines **COMPLETED** 22h11m (N=16 and N=32). Empirical Tree-KL/Split-KL = NaN. sim_neutral Tree-KL=0.693147 (ln2) for every method. TreeSBM Split-KL N=16 = **57.52**, N=32 = **122.38**; Branch W1 N=16 = 5.22e-5.

**Protocol / paths**

| Field | Value |
|---|---|
| Data | `data/covid/test` (Brazil geo); enrich groups 1–20 |
| Ckpt | `checkpoints/covid_v5_mutrec/best.pt` (ep14, val=0.9931) |
| Mask / EVEscape | `results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt`, `data/covid/evescape_spike_rbd.pt` |
| Script | `scripts/slurm_table7_refproc.sh` |
| ESM-C model | `Synthyra/ESMplusplus_small` |
| ProGen2 | `PROGEN2_HOME=/vast/projects/pranam/lab/nnori/progen2` · `progen2-small` |
| Enrich out | `checkpoints/table7_<ROW>_mrs0.5.json` |
| Coverage out | `benchmarks/results/coverage_curves_covid_N16_table7_<ROW>_K100.csv` |
| Tree-KL source | `benchmarks/results/results_baselines_covid.csv` job **7575017** |

EVE Δ = `model_evescape − gt_evescape` (all scored muts). Antigenic = `lit_hotspot_mut_frac`.

## Raw enrich summaries (no per_tree)

### JTT/WAG/LG (no fitness)  enrich=`table7_jtt_nofit_mrs0.5.json`  job 7575908

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "r0_backend": "jtt",
  "fitness_beta": 0.0,
  "ablate_bridge": true,
  "mut_recovery": 0.02256654639829583,
  "site_recall": 0.24604222296368858,
  "aa_acc_given_hit": 0.098442518712929,
  "cons_retention": 0.8208814215060034,
  "plm_nll": null,
  "identity": 0.6287797888245378,
  "dist_to_root": 217.60500000000002,
  "gt_dist_to_root": 292.1166666666667,
  "site_precision": 0.24680320669417025,
  "lit_hotspot_mut_frac": 0.0066325957243170445,
  "model_evescape": -0.023351704820169916,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": 0.8962924871792809,
  "ntd_mut_frac": 0.6039647209968122,
  "rbd_mut_frac": 0.08150675657373135,
  "mut_recovery_any_descendant": 0.3971540601614706,
  "site_recall_any_descendant": 0.9865073106845583
}
```

K=100 coverage row (`coverage_curves_covid_N16_table7_jtt_nofit_K100.csv` job 7585984):

```json
{
  "method": "treesbm",
  "N": "16",
  "K": "100",
  "n_roots": "5",
  "coverage": "1.0",
  "mean_min_edit": "1.4125",
  "site_recall": "0.6264705882352941",
  "mut_recovery": "0.00078125",
  "cons_retention": "1.0",
  "aa_acc_given_hit": "1.0",
  "coverage_obs_e0": "0.6",
  "coverage_obs_e1": "0.775",
  "coverage_obs_e2": "0.775",
  "coverage_obs_e3": "0.7875",
  "coverage_obs_e5": "0.9625"
}
```

### ESM-2-650M (no fitness)  enrich=`table7_esm2_650m_nofit_mrs0.5.json`  job 7575909

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "r0_backend": "esm2_650m",
  "fitness_beta": 0.0,
  "ablate_bridge": true,
  "mut_recovery": 0.01800009307015113,
  "site_recall": 0.16402659293648386,
  "aa_acc_given_hit": 0.10825003810695195,
  "cons_retention": 0.8931005374475165,
  "plm_nll": null,
  "identity": 0.6863400274723607,
  "dist_to_root": 125.13166666666666,
  "gt_dist_to_root": 292.1166666666667,
  "site_precision": 0.24676091434339154,
  "lit_hotspot_mut_frac": 0.005861658502272593,
  "model_evescape": -0.033717519397167335,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": 0.9558275352340353,
  "ntd_mut_frac": 0.7098815194603759,
  "rbd_mut_frac": 0.0858761630231777,
  "mut_recovery_any_descendant": 0.23695923739010696,
  "site_recall_any_descendant": 0.8788412255545757
}
```

K=100 coverage row (`coverage_curves_covid_N16_table7_esm2_650m_nofit_K100.csv` job 7585985):

```json
{
  "method": "treesbm",
  "N": "16",
  "K": "100",
  "n_roots": "5",
  "coverage": "1.0",
  "mean_min_edit": "1.4125",
  "site_recall": "0.3694117647058824",
  "mut_recovery": "0.0020833333333333333",
  "cons_retention": "1.0",
  "aa_acc_given_hit": "1.0",
  "coverage_obs_e0": "0.6",
  "coverage_obs_e1": "0.775",
  "coverage_obs_e2": "0.775",
  "coverage_obs_e3": "0.7875",
  "coverage_obs_e5": "0.975"
}
```

### ESM-2-650M (fitness)  enrich=`table7_esm2_650m_fit_mrs0.5.json`  job 7575911

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "r0_backend": "esm2_650m",
  "fitness_beta": 1.0,
  "ablate_bridge": true,
  "mut_recovery": 0.0010895220806443544,
  "site_recall": 0.019554821434337455,
  "aa_acc_given_hit": 0.09829606174192915,
  "cons_retention": 0.9883866739105919,
  "plm_nll": null,
  "identity": 0.7554276004231391,
  "dist_to_root": 19.333333333333336,
  "gt_dist_to_root": 292.1166666666667,
  "site_precision": 0.2930315985963868,
  "lit_hotspot_mut_frac": 0.0018322730273458406,
  "model_evescape": -0.4081989260530951,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": 0.9578732023847863,
  "ntd_mut_frac": 0.6755147532682196,
  "rbd_mut_frac": 0.06244658496938697,
  "mut_recovery_any_descendant": 0.03280188391621623,
  "site_recall_any_descendant": 0.3277829818601636
}
```

K=100 coverage row (`coverage_curves_covid_N16_table7_esm2_650m_fit_K100.csv` job 7585986):

```json
{
  "method": "treesbm",
  "N": "16",
  "K": "100",
  "n_roots": "5",
  "coverage": "1.0",
  "mean_min_edit": "1.425",
  "site_recall": "0.031764705882352945",
  "mut_recovery": "0.0",
  "cons_retention": "1.0",
  "aa_acc_given_hit": "nan",
  "coverage_obs_e0": "0.6",
  "coverage_obs_e1": "0.775",
  "coverage_obs_e2": "0.775",
  "coverage_obs_e3": "0.7875",
  "coverage_obs_e5": "0.9625"
}
```

### ESM-C (no fitness)  enrich=`table7_esmc_nofit_mrs0.5.json`  job **7627569** COMPLETED

Prior failures: 7585981 libcudart; 7575912 fair-esm. Rewired to `Synthyra/ESMplusplus_small`.

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "r0_backend": "esmc",
  "fitness_beta": 0.0,
  "ablate_bridge": true,
  "mut_recovery": 0.034691077105935136,
  "site_recall": 0.36305825062931707,
  "aa_acc_given_hit": 0.07960599698106258,
  "cons_retention": 0.6923237317333197,
  "plm_nll": 0.4604581510376688,
  "lit_hotspot_mut_frac": 0.005537470559108886,
  "model_evescape": 0.01017311537227285,
  "gt_evescape": 0.04712984991857271
}
```

### ESM-C (fitness)  enrich=`table7_esmc_fit_mrs0.5.json`  job **7627570** COMPLETED

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "r0_backend": "esmc",
  "fitness_beta": 1.0,
  "ablate_bridge": true,
  "mut_recovery": 0.03556098176632709,
  "site_recall": 0.3915233577537777,
  "aa_acc_given_hit": 0.08529255426807124,
  "cons_retention": 0.7133513181339514,
  "plm_nll": 0.42673972099528656,
  "lit_hotspot_mut_frac": 0.00532167925173563,
  "model_evescape": 0.013962596732167428,
  "gt_evescape": 0.04712984991857271
}
```

### ProGen2 (fitness)  enrich=`table7_progen2_fit_mrs0.5.json`  job **7627918** COMPLETED (n=20; after smoke **7627917**)

Prior **7627571** COMPLETED partial **n=3** — **superseded**. Prior 7585983: 0 trees / `get_head_mask`. Rewired enijkamp progen2-small; chunked long-seq fix for L>1023. Cov@100 job **7628262** status UNKNOWN (Day 4 Betty SSH timeout; leave Cov —).

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "r0_backend": "progen2",
  "fitness_beta": 1.0,
  "ablate_bridge": true,
  "mut_recovery": 0.11263959687095067,
  "site_recall": 0.1365019261821371,
  "aa_acc_given_hit": 0.3977601252978525,
  "cons_retention": 0.9854519578401151,
  "plm_nll": 0.4546100645571679,
  "lit_hotspot_mut_frac": 0.004798741735845365,
  "model_evescape": 0.014602785717328022,
  "gt_evescape": 0.04712984991857271
}
```

### TreeSBM default (fitness β=1)  enrich=`table7_treesbm_fit_mrs0.5.json`  job **7626625** (pLM NLL re-enrich; prior 7575915)

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "r0_backend": null,
  "fitness_beta": 1.0,
  "ablate_bridge": false,
  "mut_recovery": 0.037839179190283005,
  "site_recall": 0.10465225990530094,
  "aa_acc_given_hit": 0.2896289806923721,
  "cons_retention": 0.9411562098767925,
  "plm_nll": 0.4464431313466095,
  "identity": 0.7265394747311965,
  "lit_hotspot_mut_frac": 0.006138898104157822,
  "model_evescape": -0.046491253885236596,
  "gt_evescape": 0.04712984991857271
}
```

K=100 coverage row (`coverage_curves_covid_N16_table7_treesbm_fit_K100.csv` job 7585987 (repeat of 7575916)):

```json
{
  "method": "treesbm",
  "N": "16",
  "K": "100",
  "n_roots": "5",
  "coverage": "1.0",
  "mean_min_edit": "1.2375",
  "site_recall": "0.09411764705882353",
  "mut_recovery": "0.03833333333333334",
  "cons_retention": "1.0",
  "aa_acc_given_hit": "1.0",
  "coverage_obs_e0": "0.6",
  "coverage_obs_e1": "0.775",
  "coverage_obs_e2": "0.775",
  "coverage_obs_e3": "0.8125",
  "coverage_obs_e5": "0.975"
}
```

