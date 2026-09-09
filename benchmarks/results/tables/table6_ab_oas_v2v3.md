# Table 6 — Antibody Track C / B (Rodriguez 82)

Source: `benchmarks/results/tables/table6_ab_oas_v2v3.json`
CDR masks: PCP codon coords (7769 families)

| model | n_families | coverage_at_100_e1 | coverage_at_100_e2 | coverage_at_100_e3 | coverage_at_100_e5 | shm_load_error | cdr_mut_recall | terminal_diversity_error | lineage_rf | lineage_rf_note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| treesbm_ab_oas_v2 | 82 | 0.0061 | 0.0213 | 0.0447 | 0.0868 | 0.0778 | 0.1615 | 9.9769 | NaN | topology_forced_RF_N/A |
| treesbm_ab_oas_v3 | 82 | 0.0061 | 0.0183 | 0.0417 | 0.0874 | 0.0828 | 0.1750 | 11.1325 | NaN | topology_forced_RF_N/A |
| treesbm_ab_oas_v4 | 82 | 0.0061 | 0.0183 | 0.0417 | 0.0868 | 0.0536 | 0.1579 | 4.8946 | NaN | cosine RateHeads; fair N=20 |

## Notes

- **Coverage@100**: absolute Hamming ε; pool = up to K=100 gen leaf AAs across N=20 rollouts.
- **Terminal diversity error**: |mean pairwise Hamming(gen) − mean pairwise Hamming(gt)| (scalar companion to leaf-div W1 in primary rollout table).
- **Lineage RF**: NaN — TreeSBM/Thrifty/DASM force observed topology; CoSiNE samples lack Newick.
- **treesbm** = pathogen ckpt; **treesbm_ab** = Track B `ab_dasm_v1`; **treesbm_ab_oas** = Track A `ab_oas_1m_v1`.
- OAS samples live under `results/samples/treesbm_ab_oas` (do not overwrite `treesbm` / `treesbm_ab`).

