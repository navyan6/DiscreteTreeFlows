# `data/` — formed (train/val/test) vs TreeSBM-generated

Two different kinds of trees live under `data/`. Do not mix them.

| Role | Where | Filenames | Use for |
|------|--------|-----------|---------|
| **Formed / observed** | `data/<dataset>/{train,val,test}/` | `group_*_rooted.nwk` + `group_*_anc_aa.fasta` (+ `group_*_bl.json`) | **Training & eval** (ASR phylogenies) |
| **TreeSBM-generated** | `data/generated/{covid,h3n2,h1n1}/` | `group_*_generated*.{nwk,fasta}` | Model outputs — **not** train splits |

## Formed splits (train / val / test)

Each dataset directory has explicit split folders:

```
data/<dataset>/
  SPLIT_PROTOCOL.json    # how train/val/test were defined
  train/                 # training trees
  val/                   # validation trees
  test/                  # held-out test trees
```

Per group (same stem):

- `group_XXX_rooted.nwk` — rooted Newick (tip + ancestral node IDs)
- `group_XXX_anc_aa.fasta` — AA sequences for those nodes (required for training)
- `group_XXX_bl.json` — branch-length side data

**Tracked with AA FASTA (NCBI/INSDC):** `h1n1`, `h1n1_temporal`, `covid`, `covid_epidemic`, `hiv_geo`, `hiv_temporal`  
**Topology (+ BL) only (GISAID):** `h3n2`, `h3n2_epidemic` — no committed `anc_aa`

## TreeSBM-generated (not a train/val/test split)

```
data/generated/
  covid/   # matched + screen + seed variants
  h3n2/
  h1n1/
```

These are **model samples** from TreeSBM rollouts (`*_generated*.{nwk,fasta}`).
Same Newick↔FASTA ID join rule as formed trees, but sequences are generated, not ASR.
See `generated/README.md`. (A tiny demo subset remains under `examples/generated/`.)

## Do not use for modern training

| Path | Why |
|------|-----|
| `data/train/`, `data/validate/`, `data/test/` | Legacy flat GISAID HA dump (EPI_ISL). Prefer `data/<virus>/{train,val,test}/`. |
| Continent / leaf FASTAs (`*seqs.fasta`, `*train_group_*.fasta`) | Raw inputs; gitignored. Not the formed-tree training format. |

## Quick check

```bash
# Formed training trees for one split
ls data/h1n1/train/group_*_rooted.nwk | head

# All TreeSBM-generated outputs
ls data/generated/*/*_generated*.nwk | wc -l
```
