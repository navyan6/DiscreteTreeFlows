# Ab OAS 1M TreeSBM — postprocess fix + retrain (2026-08-13)

After **7585307 FAILED** (`Found 0 complete groups`): preprocess had
`group_XXX.nwk` + root-only `*_anc_aa.fasta`, no `*_rooted.nwk` / `*_bl.json`.

## Fix (Betty, no ANARCI redo)

Ran `scripts/ab_postprocess_clones_for_treesbm.py --data data/ab_clones_1m`:

| Split | n complete groups |
|------:|------------------:|
| train | 397 |
| val   | 44 |
| test  | 170 |

Per group: `group_XXX.nwk` → `group_XXX_rooted.nwk` (root=`NODE_ROOT`);
`*_anc_aa.fasta` = root + all leaf AAs; `*_bl.json` with
`numdate` = cumulative Newick path length from root (not zeros).

## Train

- Script: `scripts/slurm_ab_oas_treesbm_train.sh`
- `LAMBDA_BR=0` (Ab genetic BLs ≠ calendar time); entropy ON + CDR IMGT mask
- Job: **7596655** → `checkpoints/ab_oas_1m_v1`
- Prior failed: 7585307

## bl.json schema (sample)

```json
{
  "generated_by": "scripts/ab_postprocess_clones_for_treesbm.py",
  "note": "numdate = cumulative Newick path length from root (genetic distance, not calendar time)",
  "nodes": {
    "NODE_ROOT": {"numdate": 0.0, "branch_length": 0.0},
    "NODE_0000000": {"numdate": 0.0598, "branch_length": 0.0598},
    "<tip_id>": {"numdate": <cum_path>, "branch_length": <edge_bl>}
  }
}
```

TreeDataset verified: train 397 / val 44 / test 170 trees.
