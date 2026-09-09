# Antibody benchmark — full run

**Submitted (UTC):** 2026-08-11T00:13:43Z  
**Scale:** 82 families × N=20 (`configs/full.yaml`)  
**Models:** thrifty, dasm_thrifty, treesbm, cosine (+ evaluate dependency)

## Job IDs

| Model | JobID | Partition | Timelimit |
|---|---|---|---|
| thrifty | 7517208 | genoa-std-mem | 24h |
| dasm_thrifty | 7517209 | genoa-std-mem | 24h |
| treesbm | 7517210 | genoa-std-mem | 48h |
| cosine | 7517211 | b200-mig45 / mig-max | 48h |
| evaluate | 7517212 | genoa-std-mem | 4h (afterok:7517208:7517209:7517210:7517211) |

## Paths

- Frozen trees: `antibody_benchmark/data/processed/benchmark_trees.jsonl` (82)
- Samples: `antibody_benchmark/results/samples/<model>/`
- Summary: `antibody_benchmark/results/summary/`
- Logs: `/vast/projects/pranam/lab/nnori/antibody_benchmark/logs/`
- IDs JSON: `/vast/projects/pranam/lab/nnori/antibody_benchmark/full_submit_ids.json`

## Checkpoints

- TreeSBM: `~/DiscreteTreeFlows/checkpoints/best.pt`
- CoSiNE: `/vast/projects/pranam/lab/nnori/antibody_benchmark/cosine_ckpts/cosine_dasm.ckpt`

## Blockers

- Dryad deferred (optional secondary)
- CoSiNE may wait on MIG node reservation (`ReqNodeNotAvail`)
