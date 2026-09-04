# ROLE: TreeSBM-generated trees (model outputs)

**Not** train/val/test. These are TreeSBM rollout samples (Newick + AA FASTA).

Formed ASR splits for training: `data/<dataset>/{train,val,test}/` — see `../README.md`.

## Layout

```
data/generated/
  covid/   # 48 pairs
  h3n2/    # 18 pairs
  h1n1/    #  4 pairs
```

Each sample is a sibling pair:

- `group_XXX_generated*.nwk`
- `group_XXX_generated*.fasta`  (same tip/node IDs)

Filename tags:

| Tag | Meaning |
|-----|---------|
| `_matched` | Size-matched to observed leaf count (main viz pairs) |
| `_screen` | Broader screen rollouts |
| `_s42` / `_s43` / `_s44` | Extra seeds |
| (bare) `_generated` | Unmatched / alternate sample |

## Provenance

Pulled from Betty `results/{covid,h3n2,h1n1}_tree_viz{,_matched,_screen}/` (Aug 2026 viz wave).
Sequences are **model-generated**, not GISAID dumps.
