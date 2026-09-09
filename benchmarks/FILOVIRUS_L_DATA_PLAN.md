# Filovirus L — data plan (TreeSBM)

**Status:** implemented — `prepare_filo_outbreak.py` + Nextstrain outbreak ingest
**Goal:** trees that represent real transmission history; holdout that matches the 2026 BDBV forecast question.

---

## How TreeSBM uses a “tree” here

Each **group** → one FASTA + CSV → MAFFT → FastTree → augur → rooted tree → precompute → train/eval.

Each **leaf** = one genome’s L-window AA (900 aa conserved slice).  
The model learns **substitution patterns along branches** on trees whose topology comes from sequence similarity + dates.

**Therefore:** a group must be sequences that plausibly share one **outbreak / transmission chain**, not an arbitrary global date sort.

---

## Why flu-style temporal split is wrong for Ebola

| Flu / H3N2 | Filovirus L |
|------------|-------------|
| Seasonal, global circulation | **Outbreak-driven**, dead-end between epidemics |
| “Next season” = coherent antigenic drift | “Next year” may be **unrelated outbreak** on another continent |
| Many trees per year-band | One West Africa 2013–16 tree ≠ one 2018 DRC tree |
| Thousands of seqs / season | BDBV: **tens** of genomes total |

Calendar-year train/test (current `prepare_bdbv_temporal.py`) is OK for **COVID-style** “forecast 2024 from ≤2022” when N is huge. For filo it:

- Merges unrelated outbreaks into one tree (FastTree on mixed lineages → **nonsense topology**).
- Makes “test = 2020–2026 BDBV” a single 15-leaf tree while train has **one 8-leaf BDBV tree** — not enough structure for TreeSBM.

---

## Recommended grouping unit (train trees)

**Primary key:** `(filovirus_species, outbreak_id)`  
**Fallback when metadata sparse:** `(species, country, collection_year)` with max span ≤2 years and max `--group-size` 80–150.

### Outbreak labels (sources)

1. **Nextstrain ebola all-outbreaks** — `metadata.tsv` outbreak / division fields (best).
2. **Manual map** for Marburg + RESTV (NCBI collection_date + country).
3. **BDBV:** 2007 Uganda, 2012 DRC/Uganda, **2026** (Pathoplexus).

### One tree = one outbreak (usually)

| Example tree | Species | ~N (target) | Role |
|--------------|---------|-------------|------|
| WA_2013_2016 | EBOV | 1000+ (after bulk ingest) | Train |
| NordKivu_2018_2020 | EBOV | 200+ | Train or **secondary holdout** |
| SUDV_2011_Uganda | SUDV | 20–40 | Train |
| MARV_Angola_2005 | MARV | 50+ | Train |
| BDBV_2007 | BDBV | ~5 | Train (small) |
| BDBV_2012 | BDBV | ~10 | Train |
| **BDBV_2026** | BDBV | Pathoplexus | **Primary test** |

If one outbreak exceeds `--group-size`, split by **date quartile** within the same outbreak (still same epidemic), not by random accession order across outbreaks.

**Never mix species in one tree** (single-taxon strict from plan fallback ladder).

---

## Recommended splits (two tracks)

### Track A — Paper headline (sparse BDBV)

| Band | Content | Purpose |
|------|---------|---------|
| **Train** | All filovirus L except BDBV 2026: EBOV/SUDV/MARV/RESTV/TAFV outbreaks + BDBV ≤2012 | Rich RdRp dynamics |
| **Val** | One medium EBOV outbreak tree held out for early stopping (e.g. 2012 DRC or a 2021 tree) — **not empty** | Checkpoint selection |
| **Test** | **BDBV 2026 only** (+ Pathoplexus when available) | Forecast new Bundibugyo lineage |

Eval targets: L:K1738N, L:Q1770R (FJ217161 numbering).

### Track B — Method sanity (rich holdout, optional table row)

Hold out a **whole EBOV outbreak** with hundreds of genomes:

- **Test:** West Africa 2013–2016 **or** Nord-Kivu 2018–2020 (entire outbreak trees).
- **Train:** all other filo outbreaks including other EBOV eras.

This answers: “Can TreeSBM forecast variation **within a known large epidemic**?” — not the BDBV headline, but gives Coverage@K and mut_recovery with tight CIs.

---

## Species scope (pan-filo L)

| Taxon | Include | Rationale |
|-------|---------|-----------|
| EBOV, SUDV, BDBV, TAFV, RESTV | Yes | Orthoebolavirus; same L map |
| **Marburg, Ravn** | **Yes** | Same Filoviridae L architecture; well documented |
| Cuevavirus | If N≥10 | Sparse but same family |
| Mononegavirales (measles, rabies) | **No** for this paper | Different paper / pretrain-only |

---

## Data ingest (fix the ~600 EBOV cap)

Current NCBI query returns ~600 genomes; West Africa alone has **thousands** in GISAID/Nextstrain.

### Priority order

1. **Nextstrain workflow files**  
   `metadata.tsv.zst` + `alignment.fasta.zst` from `files/workflows/ebola/`  
   → accession, date, country, **outbreak**, aligned genome → slice L.

2. **NCBI broadened query** (per species + Marburg)  
   Drop `"complete genome"[Title]` requirement; keep length ≥18 kb; dedupe on accession.

3. **Pathoplexus** — BDBV 2026 test only (restricted terms).

4. **Optional:** GISAID export if Nextstrain alignment unavailable (manual).

Target after re-ingest: **≥2,000 EBOV L**, **≥50 MARV L**, **≥26 BDBV L**, many outbreak-labeled trees.

---

## Pipeline changes (before train)

| Script | Change |
|--------|--------|
| `download_ebolavirus_nextstrain.py` | Ingest workflow zst; parse outbreak column |
| `download_ebolavirus_ncbi.py` | Add MARV/RESTV; relax query |
| **`prepare_filo_outbreak.py`** (new) | Group by outbreak_id; write `SPLIT_PROTOCOL.json` with outbreak histogram |
| `prepare_bdbv_temporal.py` | Deprecate for primary track; keep as ablation |
| `slurm_bdbv_train_mut_recovery.sh` | Val = held-out outbreak dir; fail if val empty |

Output dirs:

```
data/filo_l/
  train/filo_train_group_NNN.fasta   # one outbreak per group
  val/filo_val_group_NNN.fasta       # one held-out EBOV outbreak
  test/filo_test_group_NNN.fasta     # BDBV 2026 only
  SPLIT_PROTOCOL.json
```

Checkpoints: `checkpoints/filo_l_v1_mutrec`, `filo_l_v1_lit_mutrec` (rename from bdbv_pan when stable).

---

## What each split “means” for eval

| Metric | Interpretation on filo |
|--------|-------------------------|
| **mut_recovery / site_recall** | On test outbreak tree(s), do generated leaves recover unseen substitutions along branches? |
| **Coverage@K** | Hamming coverage of generated vs true leaves in held-out outbreak |
| **K1738N / Q1770R** | Site-level forecast on BDBV 2026 test only |
| **Pan vs BDBV-only train** | Does EBOV/MARV outbreak history help BDBV 2026? |

---

## Decision checklist (before sbatch)

- [ ] Outbreak metadata joined for ≥90% of train seqs  
- [ ] ≥10 train trees, ≥50 leaves/train tree median (EBOV outbreaks)  
- [ ] Val = exactly 1–2 held-out outbreak trees (not empty)  
- [ ] Test = BDBV 2026 only (no EBOV in test for headline track)  
- [ ] No cross-species groups  
- [ ] `DATA_AUDIT.md` shows N by outbreak, not just by year  
- [ ] Pathoplexus placeholder for 2026 test refresh  

---

## Ablation matrix (unchanged intent)

| Train | Test | Question |
|-------|------|----------|
| Pan-filo outbreaks | BDBV 2026 | Primary paper |
| BDBV outbreaks only | BDBV 2026 | Does cross-species help? |
| Pan-filo | EBOV WA 2013–16 holdout | Method sanity (rich) |
| Pan-filo + lit mask | BDBV 2026 | Hotspot ablation |

Do **not** overwrite COVID/flu/HIV ckpts.
