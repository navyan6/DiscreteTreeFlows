# ROLE: TreeSBM-generated (NOT train/val/test)

These files are **sampled by TreeSBM**, not ASR-formed phylogenies.

| | Formed (train/val/test) | This folder |
|--|-------------------------|-------------|
| Path | `data/<dataset>/{train,val,test}/` | `data/examples/generated/` |
| Newick | `group_*_rooted.nwk` | `group_*_generated.nwk` |
| FASTA | `group_*_anc_aa.fasta` | `group_*_generated.fasta` |
| Sequences | Observed / ancestral reconstruction | Model-generated AA |
| Use | Training & metrics | Demos / viz only |

## Layout

```
covid/group_XXX_generated.{nwk,fasta}
h3n2/group_XXX_generated.{nwk,fasta}
```

Newick has no sequences — join tip/node IDs to the sibling FASTA.

## Provenance

Copied from `results/covid_tree_viz/` and `results/h3n2_tree_viz/`
`*_generated_matched.{nwk,fasta}` (Aug 2026 viz wave).
