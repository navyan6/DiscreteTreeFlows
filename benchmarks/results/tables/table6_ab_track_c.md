# Table 6 / Paper Table 5 — Antibody Track C (Rodriguez 82)

**Updated:** 2026-08-26. Canonical paste + experiment writeup: [`table5_ab_affinity_maturation.md`](table5_ab_affinity_maturation.md).  
Source: `t5_summaries_betty.json` + `table6_ab_track_c.json`. CDR masks: PCP codon coords. Do not invent numbers.

## Paper Table 5 paste (4 columns; Cov@100 = absolute Hamming **ε=5**, K=100)

| Method                       | CDR mut. recall ↑ | SHM load error ↓ | Terminal diversity error ↓ | Coverage@100 ↑ |
| ---------------------------- | ----------------- | ---------------- | -------------------------- | -------------- |
| **Neutral SHM** (JC69 NT)    | 0.688             | 0.054            | 17.435                     | 0.059          |
| **Autoregressive tree-edit** | 0.000             | 0.117            | 19.142                     | 0.081          |
| **CTMC (CoSiNE)**            | 0.658             | 0.037            | 3.019                      | 0.076          |
| **TreeSBM** (pathogen ckpt)  | 0.278             | 0.055            | 5.104                      | 0.084          |

### Extended paste (+ pLM)

| Method                       | CDR mut. recall ↑ | SHM load error ↓ | Terminal diversity error ↓ | Coverage@100 ↑ |
| ---------------------------- | ----------------- | ---------------- | -------------------------- | -------------- |
| **Neutral SHM** (JC69 NT)    | 0.688             | 0.054            | 17.435                     | 0.059          |
| **pLM mutation prior only**  | 0.297             | 0.032            | 2.63                       | 0.076          |
| **Autoregressive tree-edit** | 0.000             | 0.117            | 19.142                     | 0.081          |
| **CTMC (CoSiNE)**            | 0.658             | 0.037            | 3.019                      | 0.076          |
| **TreeSBM** (pathogen ckpt)  | 0.278             | 0.055            | 5.104                      | 0.084          |




### Caps / blockers


| Row              | Definition / artifact                                                                                                                                         | Status                                                                                                                                                            |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Neutral SHM**  | **JC69** independent-site NT CTMC on **observed topo+BL** (not Thrifty, not S5F neural). Samples: `results/samples/neutral_shm` (82×20).                      | **FILLED** (nohup on login01; sbatch DNS broken)                                                                                                                  |
| **pLM prior**    | ESM-2 8M on observed topo+BL (`models/plm_prior_ab.py`).                                                                                                      | **FILLED** 2026-08-16: 82×20 Rod.82 samples → `table6_ab_plm_prior_rod82.json` (Cov@e2=0.018, CDR=0.297, Cov@e5=0.076). |
| **AR tree-edit** | ARTreeFormer **N16 pool pruned** to family `n_leaves` + shared **JTT** (not Ab-trained ARTreeFormer). Free topology. Samples: `results/samples/ar_tree_edit`. | **FILLED** (tested). CDR recall **0.0** = real outcome under JTT+adapted pool, not a dash.                                                                        |
| **CoSiNE**       | Job **7517211**; unguided Gillespie on observed edges.                                                                                                        | **FILLED** (all 4 cols; was already in Track C eval)                                                                                                              |
| **TreeSBM**      | Paper Track C = **pathogen** `checkpoints/best.pt` (job **7062213** / rollout **7517210**), Rod.82 protocol.                                                  | **FILLED**. OAS `treesbm_ab_oas` (**7627919**) additive below — not the main T5 row.                                                                              |


**Neutral ≠ Thrifty.** Thrifty (`ThriftyHumV0.2-59`, job **7517208**) is kept as an additive context-SHM row only — do **not** label it Neutral SHM.

---



## Full additive Track C / B table (ε grid + extras)


| model                           | n_families | Cov@e1 | Cov@e2 | Cov@e3 | Cov@e5 | shm_load_error | cdr_mut_recall | terminal_diversity_error | lineage_rf_note                     |
| ------------------------------- | ---------- | ------ | ------ | ------ | ------ | -------------- | -------------- | ------------------------ | ----------------------------------- |
| **neutral_shm** (JC69)          | 82         | 0.0061 | 0.0152 | 0.0295 | 0.0589 | 0.0537         | 0.6883         | 17.4352                  | topology_forced                     |
| thrifty (additive; not Neutral) | 82         | 0.0061 | 0.0183 | 0.0386 | 0.0579 | 0.0366         | 0.7463         | 12.9655                  | topology_forced                     |
| dasm_thrifty                    | 82         | 0.0061 | 0.0213 | 0.0467 | 0.0854 | 0.0276         | 0.7020         | 2.7075                   | topology_forced                     |
| **ar_tree_edit**                | 82         | 0.0061 | 0.0183 | 0.0417 | 0.0813 | 0.1168         | 0.0000         | 19.1422                  | free_topo_no_Newick_RF              |
| **cosine**                      | 82         | 0.0061 | 0.0213 | 0.0447 | 0.0762 | 0.0371         | 0.678          | 3.0190                   | no_gen_newick                       |
| **treesbm** pathogen            | 82         | 0.0061 | 0.0183 | 0.0417 | 0.0843 | 0.0548         | 0.2776         | 5.1041                   | topology_forced                     |
| treesbm_ab (`ab_dasm_v1`)       | 82         | 0.0061 | 0.0213 | 0.0447 | 0.0813 | 0.0906         | 0.1157         | 12.9461                  | topology_forced                     |
| treesbm_ab_oas (`ab_oas_1m_v1`) | 82         | 0.0061 | 0.0244 | 0.0447 | 0.0813 | 0.0774         | 0.1496         | 9.8881                   | topology_forced                     |
| **plm_prior**                   | 82         | 0.0061 | 0.0183 | 0.0386 | 0.0762 | 0.0321         | 0.2969         | 2.6284                   | topology_forced                     |




## Notes

- **Coverage@100 (paper):** absolute Hamming **ε=5**; pool = up to K=100 gen leaf AAs across N=20 rollouts (same convention as prior Cov 0.058/0.084).
- **Terminal diversity error:** |mean pairwise Hamming(gen) − mean pairwise Hamming(gt)|.
- **Lineage RF:** NaN when topology forced or no gen Newick (AR free-topo Newick not exported yet).
- **treesbm** = pathogen ckpt (paper T5); **treesbm_ab_oas** = Track A OAS ckpt eval on Rod.82 (**7627919**).
- **Cluster:** BatchMode SSH → login01 where `sbatch` fails (DNS SRV `_slurmctld`); Neutral/AR/pLM launched via **nohup** on login01. Prefer `bash antibody_benchmark/scripts/betty_submit_t5_baselines.sh` from Duo **login03** when Slurm works.

