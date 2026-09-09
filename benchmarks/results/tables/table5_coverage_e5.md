# Table 5 — Coverage @ ε=5 absolute Hamming (paste-ready)

**Source:** existing e-abs CSVs (jobs **7575015** COVID / **7575016** H1). Column `coverage_obs_e5`.  
**Protocol:** N=16, K∈{10,100}, 5 diverse Brazil roots (groups 4/5/7/9/10) / H1 geo roots; ckpt `covid_v5_mutrec` / `h1n1_v2_lit_hotspot`.

**Note:** mean_min_edit ≈ 1.25–1.46 on these roots, so **ε=5 saturates** (~0.96) and methods barely separate. Prefer **ε=2** (`coverage_obs_e2`) for discrimination; ε=5 matches Table 4’s absolute Hamming definition but is too loose here.

| Virus | Method | Cov@10 | Cov@100 |
|---|---|---:|---:|
| SARS-CoV-2 Spike | pLM | 0.9625 | 0.9625 |
| SARS-CoV-2 Spike | AR | 0.9625 | 0.9625 |
| SARS-CoV-2 Spike | TreeSBM | 0.9625 | 0.9625 |
| Influenza H1N1 HA | pLM | 0.9625 | 0.9625 |
| Influenza H1N1 HA | AR | 0.9625 | 0.9625 |
| Influenza H1N1 HA | TreeSBM | **0.9750** | 0.9625 |

Companion ε=2 (same CSVs): COVID pLM/AR 0.775 / TreeSBM 0.775→0.788; H1 pLM/AR 0.762 / TreeSBM 0.800→0.762.
