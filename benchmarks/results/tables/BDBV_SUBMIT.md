# BDBV L TreeSBM — Betty submission runbook

One-page checklist when SSH to Betty is available. All scripts live under `scripts/`.

## 1. Sync repo (from laptop)

```bash
cd ~/Documents/GitHub/DiscreteTreeFlows
rsync -avz scripts/bdbv_* scripts/betty_*bdbv* scripts/download_ebolavirus_* \
  scripts/prepare_bdbv_* scripts/build_bdbv_* scripts/eval_bdbv_* \
  scripts/slurm_bdbv_* scripts/slurm_eval_bdbv.sh scripts/slurm_eval_fixed_topo.sh \
  benchmarks/BDBV_SPLITS.md benchmarks/results/tables/bdbv_* \
  nnori@login.betty.parcc.upenn.edu:~/DiscreteTreeFlows/
# Move any flattened files into scripts/ and benchmarks/ if rsync dropped them at repo root
```

Or `git pull` on Betty if changes are pushed.

## 2. Data pipeline (login node, ~2 h with NCBI refetch)

```bash
cd ~/DiscreteTreeFlows
bash scripts/betty_resume_bdbv_ingest.sh   # if ingest partially done
# OR full ingest:
bash scripts/betty_submit_bdbv.sh --ingest-only
```

Expected after ingest:

| Artifact | Notes |
|----------|-------|
| `data/bdbv/manifest.json` | ~690 NCBI genomes |
| `data/bdbv/l_nt/all.fasta` | L CDS via GenBank `gene=L` |
| `results/bdbv_l_conservation/window_config.json` | MAX_SEQ_LEN ≈ 900 |
| `data/bdbv_temporal/` | BDBV-only temporal split |
| `data/bdbv_pan_temporal/` | Pan-ebolavirus train + BDBV test |

Preflight:

```bash
bash scripts/bdbv_preflight.sh
```

## 3. Submit SLURM wave (default: no re-ingest)

```bash
bash scripts/betty_submit_bdbv.sh
```

Writes job IDs to `benchmarks/results/tables/bdbv_job_ids.json`.

### Job DAG (summary)

1. `slurm_bdbv_pipeline.sh` ×2 (bdbv_temporal + bdbv_pan_temporal)
2. `slurm_bdbv_precompute.sh` ×2 (after pipeline)
3. Train ×4: `bdbv_v1_mutrec`, `bdbv_v1_lit_mutrec`, `bdbv_pan_v1_mutrec`, `bdbv_pan_v1_lit_mutrec`
4. Eval ×4 + coverage + fixed-topo for primary checkpoint

Monitor:

```bash
squeue -u nnori
tail -f logs/bdbv_*.log
```

## 4. Pathoplexus (when access granted)

```bash
PATHOPLEXUS_SRC=/path/to/export bash scripts/bdbv_sync_pathoplexus.sh
python scripts/prepare_bdbv_temporal.py   # refresh test band
bash scripts/betty_submit_bdbv.sh --eval-only
```

## 5. Results

Fill `benchmarks/results/tables/bdbv_v1_results.md` from:

- `benchmarks/results/coverage_curves_bdbv_N16_eabs.csv`
- `checkpoints/eval_enrichment_bdbv_*_mrs*.json`
- `checkpoints/eval_bdbv_2026_variants_*.json`

Split policy: `benchmarks/BDBV_SPLITS.md`

## Known data notes (Aug 2026)

- NCBI L extract yields ~622 L CDS (GenBank `gene=L`); 2026 BDBV genomes in manifest may lack annotated L → **Pathoplexus required for true 2026 test**.
- Temporal cutoffs tuned to DATA_AUDIT: train BDBV ≤2018, test 2020–2026 (proxy until Pathoplexus).
- Nextstrain workflow URLs may 404; NCBI bulk is primary.
