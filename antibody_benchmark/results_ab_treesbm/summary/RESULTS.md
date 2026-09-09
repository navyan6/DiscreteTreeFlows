# Results summary (TEMP until full model rollouts)

## Primary (distributional rollout)

| Model | Root→leaf W1 ↓ | Site freq ρ ↑ | Substitution JS ↓ | Co-mutation ρ ↑ | Leaf diversity W1 ↓ |
| --- | --- | --- | --- | --- | --- |
| treesbm | 11.080 [10.131,12.088] | 0.144 [0.123,0.164] | 0.508 [0.494,0.520] | 0.043 [0.004,0.083] | 14.269 [13.240,15.429] |
| thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| dasm_thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |
| cosine TEMP/stub | TEMP/— | TEMP/— | TEMP/— | TEMP/— | TEMP/— |


## Secondary

| Model | Transition NLL ↓ | Branch mutation calibration ↓ | Exact descendant distance ↓ |
| --- | --- | --- | --- |
| treesbm | N/A | TEMP/— | 0.121 |
| thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| dasm_thrifty TEMP/stub | TEMP/— | TEMP/— | TEMP/— |
| cosine TEMP/stub | TEMP/— | TEMP/— | TEMP/— |


Notes:
- Recursive root-conditioned rollouts; generated parents feed children.
- Family-level bootstrap CIs when computed.
- Rows labeled TEMP/stub are null models or missing adapters.
- CoSiNE must use unguided Gillespie; TreeSBM forced to observed topology+BL.
