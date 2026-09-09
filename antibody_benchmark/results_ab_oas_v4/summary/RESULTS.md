# Results summary (TEMP until full model rollouts)

## Primary (distributional rollout)

| Model | Root→leaf W1 ↓ | Site freq ρ ↑ | Substitution JS ↓ | Co-mutation ρ ↑ | Leaf diversity W1 ↓ |
| --- | --- | --- | --- | --- | --- |
| treesbm_ab_oas_v4 | 6.624 [6.017,7.283] | -0.057 [-0.081,-0.034] | 0.464 [0.452,0.477] | 0.026 [0.009,0.043] | 7.538 [6.945,8.181] |
| thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| dasm_thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| cosine TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| treesbm TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |


## Secondary

| Model | Transition NLL ↓ | Branch mutation calibration ↓ | Exact descendant distance ↓ |
| --- | --- | --- | --- |
| treesbm_ab_oas_v4 | TEMP/— | TEMP/— | 0.148 |
| thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| dasm_thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| cosine TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| treesbm TEMP/stub | N/A | TEMP/— | TEMP/— |


Notes:
- Recursive root-conditioned rollouts; generated parents feed children.
- Family-level bootstrap CIs when computed.
- Rows labeled TEMP/stub are null models or missing adapters.
- CoSiNE must use unguided Gillespie; TreeSBM forced to observed topology+BL.
