# Filovirus L — limited data workarounds

Current Betty ingest (~1,850 windowed L seqs; **test empty** until Pathoplexus 2026) is enough to **train** TreeSBM but not for tight BDBV forecast CIs. Options ranked by paper fit.

---

## A) Data expansion (do first)

| Action | Expected gain |
|--------|----------------|
| **Nextstrain Charon ingest** (`download_ebolavirus_nextstrain.py`) | +2–4k accessions with real outbreak/division labels; sequences via NCBI efetch in `bdbv_extract_l.py` |
| **Smaller group-size (80) / min-group (5)** | More train trees from same L pool |
| **Marburg NCBI** (broader query, no title filter) | +50–200 MARV L |
| **Pathoplexus 2026 BDBV** | Gold test band for K1738N/Q1770R |
| **GISAID EBOV export** (manual) | Thousands if Charon+NCBI still thin |

---

## B) Pan-filo train (recommended now)

**What:** Train on all Filoviridae L outbreaks (EBOV + SUDV + BDBV + MARV + RESTV + TAFV) in `data/filo_l/train`; test on BDBV 2026 when available.

**Why it helps:** TreeSBM needs **branch diversity**; EBOV Kivu/WA trees teach RdRp dynamics that BDBV-only (19 seqs) cannot.

**Already implemented:** `prepare_filo_outbreak.py` + pan manifest.

**Ablation row:** BDBV-only train vs pan-filo (same test).

---

## C) Transfer learning / init (no new architecture required)

| Method | How | When to use |
|--------|-----|-------------|
| **Init from COVID/flu ckpt** | `--init-checkpoint checkpoints/covid_v5_mutrec/best.pt` + `filo_l_v1_tl_mutrec` (wired in `train.py` / submit script) | If filo train is unstable or val loss flat |
| **Two-stage** | (1) Pretrain on pan-filo + all EBOV outbreaks; (2) finetune with λ_mut↑ on BDBV historical outbreaks only | If pan train dilutes BDBV sites |
| **PLM prior only** | Already in pipeline via `precompute_plm.py`; ensure ESM embeddings cover rare BDBV columns | Default |

Implementation sketch (stage 2):
```bash
# after filo_l_v1_mutrec converges
CKPT_DIR=checkpoints/filo_l_v1_bdbv_ft \
  INIT_CKPT=checkpoints/filo_l_v1_mutrec/best.pt \
  sbatch ...  # add --init-checkpoint to train.py if wired
```

---

## D) Mononegavirales L1 pretrain (exploratory — not main paper)

**Idea:** Pretrain on **all Mononegavirales L** (measles, rabies, NiV, etc.) — thousands of seqs — then finetune filo L.

**Pros:** Maximum L1 diversity; learns generic RdRp substitution grammar.

**Cons:**
- **Site mapping breaks** for BDBV K1738N/Q1770R (different gene length / MSA gap structure).
- Different paper framing (cross-order transfer vs filovirus outbreak forecast).
- Needs new MSA window + mask rebuild.

**Verdict:** Appendix / future work unless main filo track fails. Prefer **pan-filo** over pan-L1.

---

## E) Eval workarounds while test is empty

| Track | Test | Purpose |
|-------|------|---------|
| **Track A** | BDBV 2026 (Pathoplexus) | Paper headline |
| **Track B** | `data/filo_l_track_b` — hold out WA 2013–16 or Kivu | Coverage@K + mut_recovery with N>500 |
| **Proxy** | BDBV 2020 outbreak in train; eval variant recovery on held-out leaves | Sanity until 2026 |

---

## F) Training knobs for small-N filo

- **Smaller `--group-size`** (60–80) → more trees, less mixed lineages per tree.
- **λ_mut=12, λ_cons=0.5** (covid_v5 recipe) — already in `slurm_epidemic_train_mut_recovery.sh`.
- **Lit mask** (`filo_l_v1_lit_mutrec`) — up-weight RdRp hotspots including K1738/Q1770 window.
- **Val = held-out EBOV outbreak** (auto in outbreak split) — never empty val again.

---

## Decision tree

1. Pathoplexus available? → merge test only, re-eval.
2. Pan-filo + Track B sanity passes? → run Track A BDBV eval.
3. BDBV 2026 recall still poor? → two-stage BDBV finetune (C) + lit mask ablation.
4. Still poor? → Nextstrain bulk ingest (A) before considering Mononegavirales (D).
