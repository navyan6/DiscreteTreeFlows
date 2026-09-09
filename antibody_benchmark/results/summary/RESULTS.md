# Results summary (TEMP until full model rollouts)

## Primary (distributional rollout)

| Model | Root→leaf W1 ↓ | Site freq ρ ↑ | Substitution JS ↓ | Co-mutation ρ ↑ | Leaf diversity W1 ↓ |
| --- | --- | --- | --- | --- | --- |
| cosine | 4.686 [4.189,5.202] | 0.419 [0.397,0.443] | 0.277 [0.260,0.294] | 0.006 [-0.009,0.025] | 5.523 [4.987,6.127] |
| dasm_thrifty | 3.799 [3.410,4.225] | 0.444 [0.420,0.469] | 0.277 [0.260,0.294] | 0.011 [-0.004,0.026] | 4.419 [3.977,4.952] |
| identity_null TEMP/stub | 12.476 [8.311,17.144] | 0.575 [0.458,0.684] | 0.347 [0.347,0.347] | TEMP/— | 16.933 [11.649,22.790] |
| poisson_null TEMP/stub | 4.118 [2.237,6.000] | 0.226 [0.106,0.406] | 0.563 [0.492,0.636] | -0.018 [-0.049,0.015] | 3.449 [2.028,5.134] |
| thrifty | 4.969 [4.409,5.586] | 0.284 [0.257,0.309] | 0.319 [0.303,0.334] | 0.016 [-0.001,0.032] | 8.165 [7.266,9.135] |
| treesbm | 6.746 [6.107,7.432] | 0.100 [0.082,0.117] | 0.431 [0.417,0.444] | 0.040 [0.018,0.061] | 7.752 [7.090,8.479] |


## Secondary

| Model | Transition NLL ↓ | Branch mutation calibration ↓ | Exact descendant distance ↓ |
| --- | --- | --- | --- |
| cosine | TEMP/— | TEMP/— | 0.148 |
| dasm_thrifty | TEMP/— | TEMP/— | 0.154 |
| identity_null TEMP/stub | TEMP/— | TEMP/— | 0.095 |
| poisson_null TEMP/stub | TEMP/— | TEMP/— | 0.147 |
| thrifty | TEMP/— | TEMP/— | 0.206 |
| treesbm | N/A | TEMP/— | 0.146 |


Notes:
- Recursive root-conditioned rollouts; generated parents feed children.
- Family-level bootstrap CIs when computed.
- Rows labeled TEMP/stub are null models or missing adapters.
- CoSiNE must use unguided Gillespie; TreeSBM forced to observed topology+BL.
