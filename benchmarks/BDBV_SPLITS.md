# BDBV / filovirus L-protein splits

Coordinate frame and outbreak targets for the TreeSBM Bundibugyo forecasting track.

## References

| Accession | Role |
|-----------|------|
| **NC_014373.1** | BDBV reference genome; L CDS nt 11566–18198 (2211 aa) |
| **FJ217161.1** | 2026 PHO variant numbering: **L:K1738N**, **L:Q1770R** |
| **KM034562.1** | EBOV Makona; literature **L:D759G**, GDNQ core 741–744 |

All eval sites are mapped to the training window via pan-L MSA (`results/bdbv_l_conservation/l_msa.fasta`).

## Data sources

| Source | Path | Notes |
|--------|------|-------|
| NCBI | `data/bdbv/raw/ncbi/` | Primary bulk; prefer on dedupe |
| Nextstrain | `data/bdbv/raw/nextstrain/` | Metadata enrichment |
| Pathoplexus | `data/bdbv/raw/pathoplexus_2026/` | 2026 gold test; merge-only |

Manifest: `data/bdbv/manifest.json`  
Audit: `data/bdbv/DATA_AUDIT.md`  
Split audit: `data/bdbv/FILO_SPLIT_AUDIT.json`

## Outbreak split (**primary**)

Script: `scripts/prepare_filo_outbreak.py`  
Output: `data/filo_l/{train,val,test}/`

**One tree = one outbreak** (single species). Large outbreaks split into date-ordered sub-groups of ≤120 seqs within the same `outbreak_id`.

| Band | Content |
|------|---------|
| **train** | All filovirus outbreaks except val + test (EBOV/SUDV/BDBV/MARV/RESTV/TAFV) |
| **val** | Auto: one held-out EBOV outbreak (~60 seqs) for early stopping |
| **test** | **BDBV 2026** (`bdbv_2026` outbreak) — Pathoplexus when available |

Group prefixes: `filo_train`, `filo_val`, `filo_test`

Track B (rich sanity): `--track-b --test-outbreak ebov_wa_2013_2016` → `data/filo_l_track_b/`

```bash
python scripts/prepare_filo_outbreak.py --out-base data/filo_l
python scripts/audit_filo_splits.py
```

Checkpoints: `checkpoints/filo_l_v1_mutrec`, `filo_l_v1_lit_mutrec`

See [`FILOVIRUS_L_DATA_PLAN.md`](FILOVIRUS_L_DATA_PLAN.md), [`EPIDEMIC_TREE_SPLITS.md`](EPIDEMIC_TREE_SPLITS.md).

## Temporal split (legacy ablation)

Script: `scripts/prepare_bdbv_temporal.py`  
Output: `data/bdbv_temporal/{train,val,test}/`

| Band | BDBV-only | Pan-ebolavirus (`--pan-ebolavirus`) |
|------|-----------|-------------------------------------|
| **train** | BDBV year ≤ **2018** | BDBV ≤ 2018 + EBOV/SUDV/TAFV (excl. test years) |
| **val** | BDBV **2015** (optional) | same |
| **test** | BDBV **2020–2026** | BDBV 2020–2026 only |

Group prefixes: `bdbvttrain`, `bdbvtval`, `bdbvttest`  
Default group size: 120 (tune down if BDBV N is tiny)

Pan train: `python scripts/prepare_bdbv_temporal.py --pan-ebolavirus`  
→ `data/bdbv_pan_temporal/` (override with `--out-base`)

## Geographic split (secondary)

Script: `scripts/prepare_bdbv_geo.py`  
Output: `data/bdbv_geo/`

East Africa holdout (Uganda, Sudan, South Sudan, DRC) vs rest — location generalization only.

## Fallback ladder

Use when Coverage@K or 2026 variant recall stays poor after primary runs:

1. **Pan-ebolavirus temporal** (default ablation)
2. **BDBV-only temporal** — if pan dilutes BDBV signal
3. **Smaller groups / more trees** — if lineage mixing breaks roots
4. **Single-outbreak holdout** — train 2007, test 2012, forecast 2026
5. **EBOV-only sanity** — recover **D759G** region; confirms pipeline, not paper headline

Decision rules (see `benchmarks/results/tables/bdbv_v1_results.md`):

- Temporal vs geo: higher 2026 variant recall + Coverage@K on temporal
- Pan vs BDBV-only: pan only if it lifts BDBV 2026 metrics
- Lit mask: wins only if variant + lit-stratified mut recall up, cons ≥ 0.98

## Pathoplexus merge (test-only)

When access is granted:

```bash
PATHOPLEXUS_SRC=/path/to/export bash scripts/bdbv_sync_pathoplexus.sh
python scripts/prepare_bdbv_temporal.py   # refresh test band only
# Re-run eval only — no retrain unless leakage audit fails
bash scripts/betty_submit_bdbv.sh --eval-only
```

## Training checkpoints

| Checkpoint | Scope | Lit mask |
|------------|-------|----------|
| `checkpoints/bdbv_v1_mutrec` | BDBV-only temporal | OFF |
| `checkpoints/bdbv_v1_lit_mutrec` | BDBV-only temporal | ON |
| `checkpoints/bdbv_pan_v1_mutrec` | Pan temporal | OFF |
| `checkpoints/bdbv_pan_v1_lit_mutrec` | Pan temporal | ON |

Geo secondary: `bdbv_geo_v1_mutrec`, `bdbv_pan_geo_v1_mutrec` (+ lit pairs if temporal lit wins)

Do **not** overwrite COVID/flu/HIV/OAS paper checkpoints.
