# Newick + FASTA

Each tree is a sibling pair with the same stem:

- `*.nwk` — topology (node IDs only; no sequences)
- `*.fasta` — AA sequences; headers must match Newick labels

**Formed (train/val/test):** `data/<dataset>/{train,val,test}/group_*_rooted.nwk` + `group_*_anc_aa.fasta`

**Generated (TreeSBM):** `data/generated/<virus>/{matched,screen,seeds}/group_*_generated*.{nwk,fasta}`
