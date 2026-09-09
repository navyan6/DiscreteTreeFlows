# Antibody benchmark filter summary

**Source PCP:** `/Users/navyanori/Documents/GitHub/DiscreteTreeFlows/antibody_benchmark/data/raw/dasm/extracted/dasm-experiments-data/v3/rodriguez-airr-seq-race-prod-NoWinCheck_igh_pcp_2024-11-12_MASKED_NI_noN_no-naive.csv.gz`

## Attrition by stage

| Stage | Count |
|---|---:|
| raw_families | 7769 |
| after_root_resolution | 7191 |
| after_topology_validation | 7191 |
| after_min_leaves | 82 |
| after_min_depth | 82 |
| after_productive_root | 82 |
| final_eligible | 82 |

## Exclusion reason counts

| Reason | Count |
|---|---:|
| lt_min_leaves | 7109 |
| multi_root_no_usable_naive | 578 |

**Final eligible families:** 82

## Root policy

- Prefer explicit `parent_is_naive` germline/naive node when present.
- Published DASM Zenodo v3 PCPs are upstream `no-naive` (naive edges already removed); we do **not** strip further and do **not** fabricate naive edges.
- Multi-root forests with no usable single naive/root are **excluded**.


## Source fingerprint

- path: `/Users/navyanori/Documents/GitHub/DiscreteTreeFlows/antibody_benchmark/data/raw/dasm/extracted/dasm-experiments-data/v3/rodriguez-airr-seq-race-prod-NoWinCheck_igh_pcp_2024-11-12_MASKED_NI_noN_no-naive.csv.gz`
- size_bytes: 1770898
- md5: `05719ad05f8e1e8b3bf99a8daf9b76ff`
- AB_DASM_DIR: ``
- root_source_counts: `{'topological': 82}`
