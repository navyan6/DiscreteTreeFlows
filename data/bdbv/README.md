# Filovirus L-protein (BDBV) TreeSBM track

Bundibugyo / pan-ebolavirus outbreak forecasting on the **L (RdRp) gene**.

## Quick start (Betty)

See **[benchmarks/results/tables/BDBV_SUBMIT.md](../benchmarks/results/tables/BDBV_SUBMIT.md)** for the full runbook.

```bash
bash scripts/betty_submit_bdbv.sh --ingest-only   # once: NCBI → splits
bash scripts/bdbv_preflight.sh
bash scripts/betty_submit_bdbv.sh                 # SLURM pipeline + train + eval
```

## Layout

```
data/bdbv/              manifest, raw NCBI, l_nt, l_window
data/bdbv_temporal/     BDBV-only temporal split
data/bdbv_pan_temporal/ pan-ebolavirus train + BDBV test
results/bdbv_l_*        conservation window + lit hotspot mask
checkpoints/bdbv_*      new ckpts (never overwrite COVID/flu/HIV paper dirs)
```

Split policy: [benchmarks/BDBV_SPLITS.md](../benchmarks/BDBV_SPLITS.md)
