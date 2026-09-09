# Results summary (TEMP until full model rollouts)

## Primary (distributional rollout)

| Model | Root→leaf W1 ↓ | Site freq ρ ↑ | Substitution JS ↓ | Co-mutation ρ ↑ | Leaf diversity W1 ↓ |
| --- | --- | --- | --- | --- | --- |
| treesbm_ab_oas_v3 | 10.133 [9.288,11.044] | 0.035 [0.012,0.060] | 0.428 [0.413,0.443] | 0.036 [-0.007,0.075] | 12.819 [11.901,13.833] |
| thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| dasm_thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| cosine TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| treesbm TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |


## Secondary

| Model | Transition NLL ↓ | Branch mutation calibration ↓ | Exact descendant distance ↓ |
| --- | --- | --- | --- |
| treesbm_ab_oas_v3 | TEMP/— | TEMP/— | 0.126 |
| thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| dasm_thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| cosine TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| treesbm TEMP/stub | N/A | TEMP/— | TEMP/— |


Notes:
- Recursive root-conditioned rollouts; generated parents feed children.
- Family-level bootstrap CIs when computed.
- Rows labeled TEMP/stub are null models or missing adapters.
- CoSiNE must use unguided Gillespie; TreeSBM forced to observed topology+BL.
