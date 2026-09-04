# Generated tree examples (Newick + AA FASTA)

TreeSBM **generated** trees for a few COVID Spike and H3N2 HA test groups.
**Newick does not store sequences** — use the sibling FASTA (same node/tip IDs).

## Layout

```
covid/group_XXX_generated.{nwk,fasta}
h3n2/group_XXX_generated.{nwk,fasta}
```

## Join rule

Every tip label in the Newick appears as a FASTA record id. These FASTAs are
model outputs (leaf AA sequences from TreeSBM generation), not GISAID dumps.

## Provenance

Copied from `results/covid_tree_viz/` and `results/h3n2_tree_viz/`
`*_generated_matched.{nwk,fasta}` (Aug 2026 viz wave). Demo of the
Newick+FASTA convention — not the full training split.

For **training**, use formed trees under `data/<split>/{train,val,test}/`
(`group_*_rooted.nwk` + `group_*_anc_aa.fasta` + `group_*_bl.json`).
NCBI/INSDC splits (h1n1, covid, hiv_*) ship with AA FASTAs; GISAID H3N2
ships topology (+ `bl.json`) only.
