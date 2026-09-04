# Newick + FASTA format

Tree topology and sequences are stored as a **sibling pair** with the same stem:

| File | Contents |
|------|----------|
| `*.nwk` | Rooted Newick (tip + internal node IDs). **No sequences.** |
| `*.fasta` | Amino-acid sequences; record IDs must match Newick labels |

Join rule: every tip/node label in the Newick appears as a FASTA header id.

## Formed trees (train / val / test)

```
data/<dataset>/{train,val,test}/
  group_XXX_rooted.nwk
  group_XXX_anc_aa.fasta
  group_XXX_bl.json          # optional branch-length side data
```

## TreeSBM-generated trees

```
data/generated/{covid,h3n2,h1n1}/
  group_XXX_generated*.nwk
  group_XXX_generated*.fasta
```

Filename tags: `_matched` (leaf-count matched), `_screen` (screen rollouts), `_s42`/… (seeds).
