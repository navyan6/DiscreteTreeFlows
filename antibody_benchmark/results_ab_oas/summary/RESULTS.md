# Results summary (TEMP until full model rollouts)

## Primary (distributional rollout)

| Model | Root→leaf W1 ↓ | Site freq ρ ↑ | Substitution JS ↓ | Co-mutation ρ ↑ | Leaf diversity W1 ↓ |
| --- | --- | --- | --- | --- | --- |
| treesbm_ab_oas | 9.477 [8.681,10.337] | 0.034 [0.015,0.056] | 0.466 [0.452,0.478] | 0.037 [-0.000,0.072] | 11.821 [10.931,12.746] |
| thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| dasm_thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| cosine TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| treesbm TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |


## Secondary

| Model | Transition NLL ↓ | Branch mutation calibration ↓ | Exact descendant distance ↓ |
| --- | --- | --- | --- |
| treesbm_ab_oas | TEMP/— | TEMP/— | 0.130 |
| thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| dasm_thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| cosine TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| treesbm TEMP/stub | N/A | TEMP/— | TEMP/— |


Notes:
- Recursive root-conditioned rollouts; generated parents feed children.
- Family-level bootstrap CIs when computed.
- Rows labeled TEMP/stub are null models or missing adapters.
- CoSiNE must use unguided Gillespie; TreeSBM forced to observed topology+BL.
