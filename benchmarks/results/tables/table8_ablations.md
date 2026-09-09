# Table 8 — COVID ablations (updated 2026-08-14)

ckpt=`covid_v5_mutrec` · mrs=0.5 · n_trees=20 Brazil geo-test · L=1280.
Paste companion: [`PASTE_TABLES_2026-08-14.md`](PASTE_TABLES_2026-08-14.md).

**All T8 enrich + Cov@100 COMPLETED** (tree-ctx / seq-branch / internal Cov **7597326–28**). pLM NLL re-enrich **7597319–25** COMPLETED (`summary.plm_nll`). Same 5 roots as T7; **Cov@e2 ≈ 0.775** almost everywhere (weak discriminator).

Tree-KL/Split-KL/W1 from T8b jobs (`MODE=baselines`, treesbm-only, N=16, sim_neutral). Empirical track NaN. Tree-KL saturates at ln2 for every ablation.

### Paste-ready (paper “Table 7” screenshot + antigenic)

| Variant | Tree-KL | Cons | AA\|hit | Antigenic | Cov@100 | Mut | pLM NLL |
|---|---:|---:|---:|---:|---:|---:|---:|
| w/o bridge | 0.6931 | 0.8758 | 0.0858 | 0.0043 | 0.775 | 0.0082 | 0.4465 |
| w/o fitness (not a T8 row) | — | — | — | — | — | — | — |
| w/o seq-dependent branching | — | 0.8181 | 0.2418 | 0.0114 | 0.775 | 0.0375 | 0.4540 |
| w/o BL head | 0.6931 | 0.8390 | 0.3749 | 0.0051 | 0.775 | 0.0795 | 0.4572 |
| w/o tree-context | — | 0.8448 | 0.2635 | 0.0120 | 0.775 | 0.0841 | 0.4562 |
| w/o per-site entropy | 0.6931 | 0.8346 | 0.3121 | 0.0053 | 0.775 | 0.0845 | 0.4573 |
| full | 0.6931 | 0.8389 | 0.3753 | 0.0051 | 0.775 | 0.0796 | 0.4572 |

### pLM NLL re-enrich (COMPLETED)

| Ablation | Job | plm_nll |
|---|---:|---:|
| no_bridge | 7597319 | 0.4465 |
| no_tree_ctx | 7597320 | 0.4562 |
| no_seq_branch | 7597321 | 0.4540 |
| no_bl_eval | 7597322 | 0.4572 |
| no_internal_seqs | 7597323 | 0.4565 |
| no_lit_mask | 7597325 | 0.4572 |
| no_entropy / full | 7585969 / 7585970 | 0.4573 / 0.4572 |

### Multi-ε coverage (edit-distance thresholds) — K=10 and K=100

Source: `benchmarks/results/coverage_curves_covid_N16_table8_*_K100.csv` · N=16 · 5 roots · e-abs list `{0,1,2,3,5,8,10}`.

**Wide @ K=100** (`coverage_obs_e*`):

| Ablation | e0 | e1 | e2 | e3 | e5 | e8 | e10 | min_edit@100 | mut | site | aa\|hit | Cons | pLM NLL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full | 0.600 | 0.775 | **0.775** | 0.788 | 0.963 | 0.988 | 0.988 | 1.413 | 0.0796 | 0.2189 | 0.3753 | 0.8389 | 0.4572 |
| no_bridge | 0.600 | 0.775 | **0.775** | 0.788 | 0.963 | 0.988 | 0.988 | 1.425 | 0.0082 | 0.1337 | 0.0858 | 0.8758 | 0.4465 |
| no_entropy | 0.600 | 0.775 | **0.775** | 0.813 | 0.975 | 0.988 | 0.988 | 1.225 | 0.0845 | 0.2206 | 0.3121 | 0.8346 | 0.4573 |
| no_bl_eval | 0.613 | 0.775 | **0.775** | 0.813 | 0.975 | 0.988 | 0.988 | 1.200 | 0.0795 | 0.2130 | 0.3749 | 0.8390 | 0.4572 |
| no_lit_mask | 0.625 | 0.775 | **0.788** | 0.813 | 0.975 | 0.988 | 0.988 | 1.175 | 0.0690 | 0.2088 | 0.3406 | 0.8393 | 0.4572 |
| no_tree_ctx | 0.600 | 0.775 | **0.775** | 0.813 | 0.975 | 0.988 | 0.988 | 1.225 | 0.0841 | 0.2703 | 0.2635 | 0.8448 | 0.4562 |
| no_seq_branch | 0.600 | 0.775 | **0.775** | 0.813 | 0.975 | 0.988 | 0.988 | 1.238 | 0.0375 | 0.1563 | 0.2418 | 0.8181 | 0.4540 |
| no_internal | 0.600 | 0.775 | **0.775** | 0.788 | 0.963 | 0.988 | 0.988 | 1.413 | 0.0905 | 0.3347 | 0.2701 | 0.8341 | 0.4565 |

**Cov@10 vs Cov@100 by ε** (how edit-distance coverage moves with K):

| Ablation | K | e0 | e1 | e2 | e3 | e5 | e8 | e10 | min_edit |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full | 10 | 0.600 | 0.775 | 0.775 | 0.788 | 0.963 | 0.988 | 0.988 | 1.425 |
| full | 100 | 0.600 | 0.775 | 0.775 | 0.788 | 0.963 | 0.988 | 0.988 | 1.413 |
| no_bridge | 10 | 0.600 | 0.775 | 0.775 | 0.788 | 0.963 | 0.988 | 0.988 | 1.425 |
| no_bridge | 100 | 0.600 | 0.775 | 0.775 | 0.788 | 0.963 | 0.988 | 0.988 | 1.425 |
| no_entropy | 10 | 0.600 | 0.775 | 0.775 | 0.788 | 0.963 | 0.988 | 0.988 | 1.425 |
| no_entropy | 100 | 0.600 | 0.775 | 0.775 | 0.813 | 0.975 | 0.988 | 0.988 | 1.225 |
| no_bl_eval | 10 | 0.600 | 0.775 | 0.775 | 0.813 | 0.975 | 0.988 | 0.988 | 1.225 |
| no_bl_eval | 100 | 0.613 | 0.775 | 0.775 | 0.813 | 0.975 | 0.988 | 0.988 | 1.200 |
| no_lit_mask | 10 | 0.613 | 0.775 | 0.775 | 0.788 | 0.963 | 0.988 | 0.988 | 1.400 |
| no_lit_mask | 100 | 0.625 | 0.775 | 0.788 | 0.813 | 0.975 | 0.988 | 0.988 | 1.175 |

**Readout:** coverage is flat e1→e2 (~0.775) then jumps at e3–e5 (~0.79→0.96) and saturates by e8. K=10→100 barely moves e1/e2; small gains appear at e0/e3/e5 for no_bl / no_lit / no_entropy. Ablations mainly differ in **min_edit** and mut/aa metrics, not mid-ε coverage.

| Ablation | Enrich job | Cov job | KL job | Mut recall | site_recall | aa\|hit | Cons | Antigenic (PMC) | EVE Δ | Cov@100 e1 | **Cov@100 e2** | Cov@100 e5 | min_edit@100 | Tree-KL | Split-KL | Branch W1 | pLM NLL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| TreeSBM w/o bridge matching | 7597319 | 7585972 | 7585977 | 0.0082 | 0.1337 | 0.0858 | 0.8758 | 0.0043 | -0.0644 | 0.7750 | **0.7750** | 0.9625 | 1.4250 | 0.6931 | 57.4440 | 0.000055 | 0.4465 |
| TreeSBM w/o tree-context encoder | 7597320 | 7597326 | — | 0.0841 | 0.2703 | 0.2635 | 0.8448 | 0.0120 | -0.0519 | 0.7750 | **0.7750** | 0.9750 | 1.2250 | — | — | — | 0.4562 |
| TreeSBM w/o seq-dependent branching | 7597321 | 7597327 | — | 0.0375 | 0.1563 | 0.2418 | 0.8181 | 0.0114 | — | 0.7750 | **0.7750** | 0.9750 | 1.2375 | — | — | — | 0.4540 |
| TreeSBM w/o BL head (gen-time flag) | 7597322 | 7585974 | 7585979 | 0.0795 | 0.2130 | 0.3749 | 0.8390 | 0.0051 | -0.0478 | 0.7750 | **0.7750** | 0.9750 | 1.2000 | 0.6931 | 57.5417 | 0.000030 | 0.4572 |
| TreeSBM w/o internal-node seqs | 7597323 | 7597328 | — | 0.0905 | 0.3347 | 0.2701 | 0.8341 | 0.0057 | — | 0.7750 | **0.7750** | 0.9625 | 1.4125 | — | — | — | 0.4565 |
| TreeSBM w/o lit/PMC hotspot mask | 7597325 | 7585975 | 7585980 | 0.0690 | 0.2088 | 0.3406 | 0.8393 | — | — | 0.7750 | **0.7875** | 0.9750 | 1.1750 | 0.6931 | 57.5371 | 0.000053 | 0.4572 |
| TreeSBM w/o site entropy | 7585969 | 7585973 | 7585978 | 0.0845 | 0.2206 | 0.3121 | 0.8346 | 0.0053 | -0.0584 | 0.7750 | **0.7750** | 0.9750 | 1.2250 | 0.6931 | 57.6113 | 0.000053 | 0.4573 |
| TreeSBM full | 7585970 | 7585971 | 7585976 | 0.0796 | 0.2189 | 0.3753 | 0.8389 | 0.0051 | -0.0478 | 0.7750 | **0.7750** | 0.9625 | 1.4125 | 0.6931 | 57.4722 | 0.000053 | 0.4572 |
| TreeSBM w/o BL head (retrain λ_br=0) | 7539981 | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — |

**7539981 λ_br=0 retrain:** TIMEOUT 24h at epoch 40 (`slurmstepd` 2026-08-12T11:12:54). `checkpoints/covid_ablate_no_bl/best.pt` kept: **epoch=14, val_loss=0.9984**, `mut_hotspot_mask=None`. No afterok enrichment. Do not treat as a finished ablation row.

**7575907 / 7585975 / 7585980 no_lit_mask:** enrich metrics **identical** to full (v5 trained without PMC mask; flag only drops scoring mask). Coverage@e2 slightly higher (0.7875 vs 0.775) / min_edit 1.175 vs 1.413 — generation-only noise on the same 5 roots, not a mask effect at train time.

**7585969 entropy ablation:** mut recall **0.0845** (full 0.0796), aa|hit **0.3121** (full **0.3753**) — entropy ON helps amino-acid identity given a hit. pLM NLL indistinguishable (0.4573 vs 0.4572).

**7575017** all-method COVID baselines COMPLETED. TreeSBM sim_neutral N=16: Tree-KL=0.693147, Split-KL=57.52, W1=5.22e-5 (matches T8b_full 57.47 / 5.26e-5). N=32 TreeSBM Split-KL=122.38, W1=8.08e-5. pLM NLL not in CSV.

**Cov done for all ablations** (7597326–28). **KL still missing** for tree-ctx / seq-branch / internal-seqs (not discriminative anyway — Tree-KL=ln2).

## Raw enrich summaries

### TreeSBM w/o bridge matching  `table8_no_bridge_mrs0.5.json`  flags=`--ablate-bridge`

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "ablate_bridge": null,
  "ablate_tree_context": null,
  "ablate_branch_length_head": null,
  "ablate_internal_node_seqs": null,
  "ablate_site_entropy": null,
  "no_lit_hotspot_mask": null,
  "branching_mode": null,
  "mut_recovery": 0.008213127596658367,
  "site_recall": 0.13374774607470574,
  "aa_acc_given_hit": 0.08583798935359593,
  "cons_retention": 0.8760794741715714,
  "plm_nll": null,
  "identity": 0.671254700553574,
  "dist_to_root": 155.36333333333332,
  "site_precision": 0.2472978568276949,
  "lit_hotspot_mut_frac": 0.004328054955664797,
  "model_evescape": -0.017254381559796896,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": 0.9737116776499747,
  "ntd_mut_frac": 0.7320912683996281,
  "rbd_mut_frac": 0.05681779919986944,
  "mut_recovery_any_descendant": 0.22788203079563027
}
```

K=100 `coverage_curves_covid_N16_table8_no_bridge_K100.csv`:

```json
{
  "method": "treesbm",
  "N": "16",
  "K": "100",
  "n_roots": "5",
  "coverage": "1.0",
  "coverage_obs_e0": "0.6",
  "coverage_obs_e1": "0.775",
  "coverage_obs_e2": "0.775",
  "coverage_obs_e3": "0.7875",
  "coverage_obs_e5": "0.9625",
  "mean_min_edit": "1.425",
  "site_recall": "0.1776470588235294",
  "mut_recovery": "0.0",
  "aa_acc_given_hit": "nan"
}
```

baselines mean sim_neutral N=16 `results_baselines_covid_t8_no_bridge.csv`: `{"tree_kl": 0.6931471609599456, "split_kl": 57.443971942377075, "branch_w_all": 5.492900784976274e-05, "rf": 0.8524248120300755, "terminal_edit": 0.00038151305453342703}`

### TreeSBM w/o tree-context encoder  `table8_no_tree_ctx_mrs0.5.json`  flags=`--ablate-tree-context`

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "ablate_bridge": null,
  "ablate_tree_context": null,
  "ablate_branch_length_head": null,
  "ablate_internal_node_seqs": null,
  "ablate_site_entropy": null,
  "no_lit_hotspot_mask": null,
  "branching_mode": null,
  "mut_recovery": 0.08407148684360902,
  "site_recall": 0.27032585544110954,
  "aa_acc_given_hit": 0.2635466180976892,
  "cons_retention": 0.8447558029036827,
  "plm_nll": null,
  "identity": 0.6552493551316071,
  "dist_to_root": 219.64666666666668,
  "site_precision": 0.2659789361938504,
  "lit_hotspot_mut_frac": 0.012011685134182971,
  "model_evescape": -0.004791226140215089,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": 0.833732010783,
  "ntd_mut_frac": 0.39621524756065973,
  "rbd_mut_frac": 0.14638197626755328,
  "mut_recovery_any_descendant": 0.3881963121531169
}
```

### TreeSBM w/o seq-dependent branching  `table8_no_seq_branch_mrs0.5.json`  flags=`poisson_ref --ref-lambda 1.0`

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "ablate_bridge": null,
  "ablate_tree_context": null,
  "ablate_branch_length_head": null,
  "ablate_internal_node_seqs": null,
  "ablate_site_entropy": null,
  "no_lit_hotspot_mask": null,
  "branching_mode": null,
  "mut_recovery": 0.045173818431861856,
  "site_recall": 0.17927069123281408,
  "aa_acc_given_hit": 0.24715213537962807,
  "cons_retention": 0.8240408112837233,
  "plm_nll": null,
  "identity": 0.6365201411522253,
  "dist_to_root": 233.3383333333333,
  "site_precision": 0.2520936835334556,
  "lit_hotspot_mut_frac": 0.01170312602455864,
  "model_evescape": 0.019221047645797264,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": 0.9288020461346164,
  "ntd_mut_frac": 0.32093782183109715,
  "rbd_mut_frac": 0.1592505040171392,
  "mut_recovery_any_descendant": 0.045173818431861856
}
```

### TreeSBM w/o BL head (gen-time flag)  `table8_no_bl_eval_mrs0.5.json`  flags=`--ablate-branch-length-head`

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "ablate_bridge": null,
  "ablate_tree_context": null,
  "ablate_branch_length_head": null,
  "ablate_internal_node_seqs": null,
  "ablate_site_entropy": null,
  "no_lit_hotspot_mask": null,
  "branching_mode": null,
  "mut_recovery": 0.0794715474850733,
  "site_recall": 0.21303103934551765,
  "aa_acc_given_hit": 0.3749175076680972,
  "cons_retention": 0.8389876252632226,
  "plm_nll": null,
  "identity": 0.6563303328199264,
  "dist_to_root": 230.525,
  "site_precision": 0.2684570267744675,
  "lit_hotspot_mut_frac": 0.005144428266790912,
  "model_evescape": -0.0007088171036399325,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": 0.9099179341296577,
  "ntd_mut_frac": 0.7206762716198417,
  "rbd_mut_frac": 0.06540384547176459,
  "mut_recovery_any_descendant": 0.5708076204650131
}
```

K=100 `coverage_curves_covid_N16_table8_no_bl_eval_K100.csv`:

```json
{
  "method": "treesbm",
  "N": "16",
  "K": "100",
  "n_roots": "5",
  "coverage": "1.0",
  "coverage_obs_e0": "0.6125",
  "coverage_obs_e1": "0.775",
  "coverage_obs_e2": "0.775",
  "coverage_obs_e3": "0.8125",
  "coverage_obs_e5": "0.975",
  "mean_min_edit": "1.2",
  "site_recall": "0.20941176470588235",
  "mut_recovery": "0.06451140873015873",
  "aa_acc_given_hit": "1.0"
}
```

baselines mean sim_neutral N=16 `results_baselines_covid_t8_no_bl_eval.csv`: `{"tree_kl": 0.6931471610349456, "split_kl": 57.541662374147016, "branch_w_all": 2.954650443231443e-05, "rf": 0.8542199248120301, "terminal_edit": 0.00043754134452391786}`

### TreeSBM w/o internal-node seqs  `table8_no_internal_seqs_mrs0.5.json`  flags=`--ablate-internal-node-seqs`

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "ablate_bridge": null,
  "ablate_tree_context": null,
  "ablate_branch_length_head": null,
  "ablate_internal_node_seqs": null,
  "ablate_site_entropy": null,
  "no_lit_hotspot_mask": null,
  "branching_mode": null,
  "mut_recovery": 0.09133459389989743,
  "site_recall": 0.356280417121234,
  "aa_acc_given_hit": 0.23466830247119227,
  "cons_retention": 0.8327929787214376,
  "plm_nll": null,
  "identity": 0.6501207661870015,
  "dist_to_root": 237.61333333333332,
  "site_precision": 0.2808911264871684,
  "lit_hotspot_mut_frac": 0.005655272984125197,
  "model_evescape": -0.014309776074320418,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": 0.9230967487965789,
  "ntd_mut_frac": 0.687870656373603,
  "rbd_mut_frac": 0.07367389953721387,
  "mut_recovery_any_descendant": 0.44205827351944815
}
```

### TreeSBM w/o lit/PMC hotspot mask  `table8_no_lit_mask_mrs0.5.json`  flags=`--no-lit-hotspot-mask`

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "ablate_bridge": false,
  "ablate_tree_context": null,
  "ablate_branch_length_head": null,
  "ablate_internal_node_seqs": null,
  "ablate_site_entropy": null,
  "no_lit_hotspot_mask": true,
  "branching_mode": null,
  "mut_recovery": 0.07956162637113982,
  "site_recall": 0.21889215045662874,
  "aa_acc_given_hit": 0.37531169853946655,
  "cons_retention": 0.8388891462612198,
  "plm_nll": null,
  "identity": 0.6563068838714516,
  "dist_to_root": 230.675,
  "site_precision": 0.26856233168989896,
  "lit_hotspot_mut_frac": null,
  "model_evescape": -0.0006701228697724671,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": null,
  "ntd_mut_frac": 0.7206863337488377,
  "rbd_mut_frac": 0.0653958477141862,
  "mut_recovery_any_descendant": 0.5710478308278572
}
```

K=100 `coverage_curves_covid_N16_table8_no_lit_mask_K100.csv`:

```json
{
  "method": "treesbm",
  "N": "16",
  "K": "100",
  "n_roots": "5",
  "coverage": "1.0",
  "coverage_obs_e0": "0.625",
  "coverage_obs_e1": "0.775",
  "coverage_obs_e2": "0.7875",
  "coverage_obs_e3": "0.8125",
  "coverage_obs_e5": "0.975",
  "mean_min_edit": "1.175",
  "site_recall": "0.5211764705882354",
  "mut_recovery": "0.11191881613756614",
  "aa_acc_given_hit": "1.0"
}
```

baselines mean sim_neutral N=16 `results_baselines_covid_t8_no_lit_mask.csv`: `{"tree_kl": 0.6931471610349456, "split_kl": 57.537074193419684, "branch_w_all": 5.253045220774457e-05, "rf": 0.8594454887218046, "terminal_edit": 0.0004450867201389177}`

### TreeSBM w/o site entropy  `table8_no_entropy_mrs0.5.json`  flags=`--ablate-site-entropy`

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "ablate_bridge": false,
  "ablate_tree_context": false,
  "ablate_branch_length_head": false,
  "ablate_internal_node_seqs": false,
  "ablate_site_entropy": true,
  "no_lit_hotspot_mask": false,
  "branching_mode": "learned",
  "mut_recovery": 0.08445981936695668,
  "site_recall": 0.22061705829621608,
  "aa_acc_given_hit": 0.312101941013336,
  "cons_retention": 0.8346305373250692,
  "plm_nll": 0.4572852925872093,
  "identity": 0.6539274711134148,
  "dist_to_root": 236.39833333333326,
  "site_precision": 0.25899611628374475,
  "lit_hotspot_mut_frac": 0.005271650104308192,
  "model_evescape": -0.011318540274749635,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": 0.8990729730333111,
  "ntd_mut_frac": 0.7075552267493672,
  "rbd_mut_frac": 0.06863593547001286,
  "mut_recovery_any_descendant": 0.6022042216029148
}
```

K=100 `coverage_curves_covid_N16_table8_no_entropy_K100.csv`:

```json
{
  "method": "treesbm",
  "N": "16",
  "K": "100",
  "n_roots": "5",
  "coverage": "1.0",
  "coverage_obs_e0": "0.6",
  "coverage_obs_e1": "0.775",
  "coverage_obs_e2": "0.775",
  "coverage_obs_e3": "0.8125",
  "coverage_obs_e5": "0.975",
  "mean_min_edit": "1.225",
  "site_recall": "0.35764705882352943",
  "mut_recovery": "0.039114583333333335",
  "aa_acc_given_hit": "1.0"
}
```

baselines mean sim_neutral N=16 `results_baselines_covid_t8_no_entropy.csv`: `{"tree_kl": 0.6931471610349456, "split_kl": 57.611306542404954, "branch_w_all": 5.312302898193719e-05, "rf": 0.8490601503759398, "terminal_edit": 0.0004389625625335924}`

### TreeSBM full  `table8_full_mrs0.5.json`  flags=`(none; 7585970 re-enrich + pLM NLL)`

```json
{
  "checkpoint": "checkpoints/covid_v5_mutrec/best.pt",
  "n_trees": 20,
  "ablate_bridge": false,
  "ablate_tree_context": false,
  "ablate_branch_length_head": false,
  "ablate_internal_node_seqs": false,
  "ablate_site_entropy": false,
  "no_lit_hotspot_mask": false,
  "branching_mode": "learned",
  "mut_recovery": 0.07956162637113982,
  "site_recall": 0.21889215045662874,
  "aa_acc_given_hit": 0.37531169853946655,
  "cons_retention": 0.8388891462612198,
  "plm_nll": 0.45718034220243053,
  "identity": 0.6563068838714516,
  "dist_to_root": 230.675,
  "site_precision": 0.26856233168989896,
  "lit_hotspot_mut_frac": 0.005146331899057694,
  "model_evescape": -0.0006701228697724671,
  "gt_evescape": 0.04712984991857271,
  "evescape_mean_antigenic_muts": 0.9118525351499324,
  "ntd_mut_frac": 0.7206863337488377,
  "rbd_mut_frac": 0.0653958477141862,
  "mut_recovery_any_descendant": 0.5710478308278572
}
```

K=100 `coverage_curves_covid_N16_table8_full_K100.csv`:

```json
{
  "method": "treesbm",
  "N": "16",
  "K": "100",
  "n_roots": "5",
  "coverage": "1.0",
  "coverage_obs_e0": "0.6",
  "coverage_obs_e1": "0.775",
  "coverage_obs_e2": "0.775",
  "coverage_obs_e3": "0.7875",
  "coverage_obs_e5": "0.9625",
  "mean_min_edit": "1.4125",
  "site_recall": "0.25588235294117645",
  "mut_recovery": "0.00078125",
  "aa_acc_given_hit": "1.0"
}
```

baselines mean sim_neutral N=16 `results_baselines_covid_t8_full.csv`: `{"tree_kl": 0.6931471610349456, "split_kl": 57.47223533011551, "branch_w_all": 5.257753394145853e-05, "rf": 0.848063909774436, "terminal_edit": 0.00044049360193492375}`

