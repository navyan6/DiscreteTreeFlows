# VOC threat tree screen

| gid | dates | n_leaves | best_voc | prospective | exact_recall | bundle_frac | root_bundle |
|-----|-------|----------|----------|-------------|--------------|-------------|-------------|
| 034 | 2020-03-11…2022-06-24 | 29 | Delta | True | 0.83 | 0.000 | [] |
| 289 | 2020-03-10…2022-07-19 | 163 | Zeta | True | 1.00 | 0.012 | [] |
| 188 | 2020-01-27…2020-06-05 | 300 | Kappa | True | 0.50 | 0.000 | [] |
| 189 | 2020-06-05…2020-07-06 | 297 | Lambda | True | 0.67 | 0.000 | [] |
| 190 | 2020-07-06…2021-02-25 | 300 | Kappa | True | 1.00 | 0.013 | ['L452R'] |
| 191 | 2021-02-28…2021-06-03 | 299 | Delta | True | 1.00 | 0.381 | ['L452R'] |
| 192 | 2021-06-04…2021-12-23 | 300 | Delta | True | 1.00 | 0.717 | ['T478K'] |
| 193 | 2021-12-23…2022-12-24 | 300 | Delta | True | 1.00 | 0.017 | ['T478K'] |
| 300 | — | — | ERROR | — | — | — | missing data/covid/train/group_300_anc_aa.fasta |
| 301 | 2020-07-03…2020-10-20 | 300 | Zeta | True | 1.00 | 0.020 | [] |
| 302 | 2020-10-20…2021-01-21 | 299 | Beta | True | 1.00 | 0.064 | ['E484K'] |
| 303 | 2021-01-24…2021-06-15 | 300 | Delta | True | 1.00 | 0.037 | ['L452R'] |
| 304 | — | — | ERROR | — | — | — | missing data/covid/train/group_304_anc_aa.fasta |
| 308 | 2021-09-14…2021-12-23 | 300 | Delta | True | 1.00 | 0.117 | ['L452R'] |
| 309 | 2021-12-26…2023-01-31 | 280 | Delta | True | 0.83 | 0.000 | ['L452R'] |
| 323 | 2020-01-29…2026-03-04 | 206 | Omicron_BA1 | True | 1.00 | 0.010 | [] |

## Per-VOC highlights (prospective + recall≥0.3)

### Gamma
- group 289 (2020-03-10…2022-07-19): exact_recall=1.00 bundle_frac=0.006 acquired_recall=0.25
- group 190 (2020-07-06…2021-02-25): exact_recall=0.42 bundle_frac=0.000 acquired_recall=0.33
- group 191 (2021-02-28…2021-06-03): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.25
- group 192 (2021-06-04…2021-12-23): exact_recall=0.42 bundle_frac=0.000 acquired_recall=0.25
- group 301 (2020-07-03…2020-10-20): exact_recall=0.42 bundle_frac=0.000 acquired_recall=0.42
- group 302 (2020-10-20…2021-01-21): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.33
- group 303 (2021-01-24…2021-06-15): exact_recall=0.58 bundle_frac=0.000 acquired_recall=0.33
- group 308 (2021-09-14…2021-12-23): exact_recall=0.42 bundle_frac=0.000 acquired_recall=0.42
- group 323 (2020-01-29…2026-03-04): exact_recall=0.42 bundle_frac=0.000 acquired_recall=0.08

### Zeta
- group 034 (2020-03-11…2022-06-24): exact_recall=0.67 bundle_frac=0.034 acquired_recall=0.33
- group 289 (2020-03-10…2022-07-19): exact_recall=1.00 bundle_frac=0.012 acquired_recall=0.00
- group 188 (2020-01-27…2020-06-05): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.00
- group 189 (2020-06-05…2020-07-06): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.00
- group 190 (2020-07-06…2021-02-25): exact_recall=0.67 bundle_frac=0.000 acquired_recall=0.33
- group 191 (2021-02-28…2021-06-03): exact_recall=0.67 bundle_frac=0.007 acquired_recall=0.33
- group 192 (2021-06-04…2021-12-23): exact_recall=0.67 bundle_frac=0.007 acquired_recall=0.33
- group 193 (2021-12-23…2022-12-24): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.00
- group 301 (2020-07-03…2020-10-20): exact_recall=1.00 bundle_frac=0.020 acquired_recall=1.00
- group 303 (2021-01-24…2021-06-15): exact_recall=0.67 bundle_frac=0.033 acquired_recall=0.33
- group 308 (2021-09-14…2021-12-23): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.33
- group 323 (2020-01-29…2026-03-04): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.00

### Alpha
- group 289 (2020-03-10…2022-07-19): exact_recall=0.43 bundle_frac=0.000 acquired_recall=0.00
- group 188 (2020-01-27…2020-06-05): exact_recall=0.43 bundle_frac=0.000 acquired_recall=0.29
- group 189 (2020-06-05…2020-07-06): exact_recall=0.57 bundle_frac=0.000 acquired_recall=0.43
- group 190 (2020-07-06…2021-02-25): exact_recall=0.86 bundle_frac=0.000 acquired_recall=0.71
- group 191 (2021-02-28…2021-06-03): exact_recall=0.86 bundle_frac=0.000 acquired_recall=0.71
- group 192 (2021-06-04…2021-12-23): exact_recall=0.57 bundle_frac=0.053 acquired_recall=0.29
- group 193 (2021-12-23…2022-12-24): exact_recall=0.71 bundle_frac=0.223 acquired_recall=0.43
- group 301 (2020-07-03…2020-10-20): exact_recall=0.57 bundle_frac=0.000 acquired_recall=0.43
- group 302 (2020-10-20…2021-01-21): exact_recall=0.43 bundle_frac=0.000 acquired_recall=0.29
- group 303 (2021-01-24…2021-06-15): exact_recall=0.57 bundle_frac=0.000 acquired_recall=0.43
- group 308 (2021-09-14…2021-12-23): exact_recall=0.57 bundle_frac=0.000 acquired_recall=0.57
- group 323 (2020-01-29…2026-03-04): exact_recall=0.57 bundle_frac=0.010 acquired_recall=0.00

### Beta
- group 034 (2020-03-11…2022-06-24): exact_recall=0.43 bundle_frac=0.000 acquired_recall=0.29
- group 289 (2020-03-10…2022-07-19): exact_recall=0.57 bundle_frac=0.000 acquired_recall=0.00
- group 189 (2020-06-05…2020-07-06): exact_recall=0.57 bundle_frac=0.000 acquired_recall=0.43
- group 190 (2020-07-06…2021-02-25): exact_recall=0.71 bundle_frac=0.000 acquired_recall=0.57
- group 191 (2021-02-28…2021-06-03): exact_recall=0.71 bundle_frac=0.000 acquired_recall=0.57
- group 192 (2021-06-04…2021-12-23): exact_recall=0.71 bundle_frac=0.000 acquired_recall=0.43
- group 193 (2021-12-23…2022-12-24): exact_recall=0.86 bundle_frac=0.000 acquired_recall=0.57
- group 301 (2020-07-03…2020-10-20): exact_recall=0.86 bundle_frac=0.000 acquired_recall=0.86
- group 302 (2020-10-20…2021-01-21): exact_recall=1.00 bundle_frac=0.064 acquired_recall=0.71
- group 303 (2021-01-24…2021-06-15): exact_recall=1.00 bundle_frac=0.027 acquired_recall=0.86
- group 308 (2021-09-14…2021-12-23): exact_recall=0.86 bundle_frac=0.000 acquired_recall=0.86
- group 309 (2021-12-26…2023-01-31): exact_recall=0.57 bundle_frac=0.000 acquired_recall=0.43
- group 323 (2020-01-29…2026-03-04): exact_recall=0.43 bundle_frac=0.000 acquired_recall=0.00

### Delta
- group 034 (2020-03-11…2022-06-24): exact_recall=0.83 bundle_frac=0.000 acquired_recall=0.67
- group 289 (2020-03-10…2022-07-19): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.00
- group 188 (2020-01-27…2020-06-05): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.17
- group 189 (2020-06-05…2020-07-06): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.33
- group 190 (2020-07-06…2021-02-25): exact_recall=0.83 bundle_frac=0.000 acquired_recall=0.33
- group 191 (2021-02-28…2021-06-03): exact_recall=1.00 bundle_frac=0.381 acquired_recall=0.50
- group 192 (2021-06-04…2021-12-23): exact_recall=1.00 bundle_frac=0.717 acquired_recall=0.50
- group 193 (2021-12-23…2022-12-24): exact_recall=1.00 bundle_frac=0.017 acquired_recall=0.50
- group 301 (2020-07-03…2020-10-20): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.33
- group 302 (2020-10-20…2021-01-21): exact_recall=0.83 bundle_frac=0.000 acquired_recall=0.50
- group 303 (2021-01-24…2021-06-15): exact_recall=1.00 bundle_frac=0.037 acquired_recall=0.67
- group 308 (2021-09-14…2021-12-23): exact_recall=1.00 bundle_frac=0.117 acquired_recall=0.50
- group 309 (2021-12-26…2023-01-31): exact_recall=0.83 bundle_frac=0.000 acquired_recall=0.50
- group 323 (2020-01-29…2026-03-04): exact_recall=1.00 bundle_frac=0.005 acquired_recall=0.17

### Omicron_BA1
- group 191 (2021-02-28…2021-06-03): exact_recall=0.36 bundle_frac=0.000 acquired_recall=0.29
- group 192 (2021-06-04…2021-12-23): exact_recall=1.00 bundle_frac=0.053 acquired_recall=0.46
- group 303 (2021-01-24…2021-06-15): exact_recall=0.39 bundle_frac=0.000 acquired_recall=0.32
- group 308 (2021-09-14…2021-12-23): exact_recall=0.39 bundle_frac=0.000 acquired_recall=0.36
- group 309 (2021-12-26…2023-01-31): exact_recall=0.32 bundle_frac=0.000 acquired_recall=0.18
- group 323 (2020-01-29…2026-03-04): exact_recall=1.00 bundle_frac=0.010 acquired_recall=0.07

### Lambda
- group 034 (2020-03-11…2022-06-24): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.17
- group 289 (2020-03-10…2022-07-19): exact_recall=1.00 bundle_frac=0.006 acquired_recall=0.00
- group 189 (2020-06-05…2020-07-06): exact_recall=0.67 bundle_frac=0.000 acquired_recall=0.50
- group 190 (2020-07-06…2021-02-25): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.33
- group 191 (2021-02-28…2021-06-03): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.33
- group 192 (2021-06-04…2021-12-23): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.33
- group 193 (2021-12-23…2022-12-24): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.17
- group 301 (2020-07-03…2020-10-20): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.33
- group 303 (2021-01-24…2021-06-15): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.17
- group 308 (2021-09-14…2021-12-23): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.33
- group 323 (2020-01-29…2026-03-04): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.00

### Mu
- group 034 (2020-03-11…2022-06-24): exact_recall=0.56 bundle_frac=0.000 acquired_recall=0.44
- group 289 (2020-03-10…2022-07-19): exact_recall=0.67 bundle_frac=0.006 acquired_recall=0.00
- group 189 (2020-06-05…2020-07-06): exact_recall=0.44 bundle_frac=0.000 acquired_recall=0.33
- group 190 (2020-07-06…2021-02-25): exact_recall=0.67 bundle_frac=0.000 acquired_recall=0.56
- group 191 (2021-02-28…2021-06-03): exact_recall=0.78 bundle_frac=0.000 acquired_recall=0.56
- group 192 (2021-06-04…2021-12-23): exact_recall=0.78 bundle_frac=0.027 acquired_recall=0.33
- group 193 (2021-12-23…2022-12-24): exact_recall=0.78 bundle_frac=0.000 acquired_recall=0.44
- group 301 (2020-07-03…2020-10-20): exact_recall=0.67 bundle_frac=0.000 acquired_recall=0.67
- group 302 (2020-10-20…2021-01-21): exact_recall=0.44 bundle_frac=0.064 acquired_recall=0.11
- group 303 (2021-01-24…2021-06-15): exact_recall=0.67 bundle_frac=0.027 acquired_recall=0.56
- group 308 (2021-09-14…2021-12-23): exact_recall=0.44 bundle_frac=0.000 acquired_recall=0.44
- group 309 (2021-12-26…2023-01-31): exact_recall=0.33 bundle_frac=0.000 acquired_recall=0.11
- group 323 (2020-01-29…2026-03-04): exact_recall=0.78 bundle_frac=0.000 acquired_recall=0.11

### Kappa
- group 034 (2020-03-11…2022-06-24): exact_recall=0.75 bundle_frac=0.000 acquired_recall=0.50
- group 289 (2020-03-10…2022-07-19): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.00
- group 188 (2020-01-27…2020-06-05): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.25
- group 189 (2020-06-05…2020-07-06): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.25
- group 190 (2020-07-06…2021-02-25): exact_recall=1.00 bundle_frac=0.013 acquired_recall=0.25
- group 192 (2021-06-04…2021-12-23): exact_recall=1.00 bundle_frac=0.007 acquired_recall=0.75
- group 193 (2021-12-23…2022-12-24): exact_recall=1.00 bundle_frac=0.003 acquired_recall=0.50
- group 301 (2020-07-03…2020-10-20): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.25
- group 302 (2020-10-20…2021-01-21): exact_recall=0.75 bundle_frac=0.000 acquired_recall=0.50
- group 303 (2021-01-24…2021-06-15): exact_recall=0.75 bundle_frac=0.000 acquired_recall=0.25
- group 308 (2021-09-14…2021-12-23): exact_recall=0.75 bundle_frac=0.000 acquired_recall=0.25
- group 309 (2021-12-26…2023-01-31): exact_recall=0.50 bundle_frac=0.000 acquired_recall=0.25
- group 323 (2020-01-29…2026-03-04): exact_recall=1.00 bundle_frac=0.000 acquired_recall=0.00

