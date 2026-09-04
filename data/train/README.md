# ROLE: legacy GISAID HA dump — not the modern split layout

Prefer `data/<virus>/{train,val,test}/` (e.g. `data/h1n1/train/`).  
This flat `data/train/` folder mixes EPI_ISL groups and is not the canonical train/val/test API.
