# Nat Microbiol 2016 → TreeSBM flu HA modeling notes

Source: Li et al., *Nat Microbiol* **1**, 16058 (2016). DOI [10.1038/nmicrobiol.2016.58](https://doi.org/10.1038/nmicrobiol.2016.58) (PMID 27572841 / PMC5087998).

## Accessibility of the Cursor upload

**`uploads/nmicrobiol201658-0.md` is a Nature paywall stub**, not the full article.

Present in the upload: title/authors, abstract, “Access options / Subscribe”, figure *titles*, acknowledgements, author list, competing interests, **supplementary table titles** (not XLSX bodies), references, HTML chrome.

**Missing from the upload:** Results/Discussion body, Methods, figure legends with genotypes, and all supplementary mutation tables.

Site lists in `src/flu_lit_sites.py` were curated from the **Cambridge institutional author manuscript** of the **same DOI** (author-accepted text + Fig. 6 genotype list). Supplementary XLSX tables were **not** downloaded → frequency-complete escape catalogs remain a **GAP**. Do not treat upload-parsed regex tokens as a site list.

## STRICT SEPARATION (flu ≠ COVID)

| | Flu (this dir) | COVID |
|---|---|---|
| Artifacts | `results/flu_mutfreq_vs_lit/` | `results/covid_mutfreq_vs_lit/` |
| Lit mask | `mut_hotspot_mask_nmicrobiol_lit.pt` (alias `mut_hotspot_mask_flu_lit.pt`) | `mut_hotspot_mask_pmc_lit.pt` |
| Region band | `mut_hotspot_mask_h3_globular_head.pt` (H3 #63–252) | `mut_hotspot_mask_rbd_band.pt` |
| Coords | H3 mature numbering + 16-aa signal → full-ORF cols | Wuhan / P0DTC2 spike |
| Eval metric | `flu_hotspot_mut_frac` / `gt_flu_hotspot_mut_frac` (+ `lit_*`) | `pmc_hotspot_mut_frac` / `gt_pmc_hotspot_mut_frac` (+ `lit_*`) |
| Site module | `src/flu_lit_sites.py` only | `PMC10142771_*` in `compute_msa_mut_freq_vs_lit.py` only |

**Rules enforced:**

1. Flu train/eval **never** loads `mut_hotspot_mask_pmc_lit.pt` or spike/RBD coords.
2. COVID train/eval **never** uses Koel / H3 / Sa flu sites.
3. Shared `eval_evescape_enrichment.py` takes an **explicit** `--lit-hotspot-mask` path with **no length-based cross-default**. Cross-pathogen path+L mismatches raise `SystemExit`.
4. COVID slurm scripts pass `results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt` only.

## Numbering (H3)

- Paper / WHO: **H3 numbering** (mature HA1; residue 1 after signal).
- Our TreeSBM HA: **full-length ungapped ORF** `L=566` including 16-aa signal `MKTIIALSYILCLVFA`, mature start `QKLPGNDNSTATLCL…`.
- Map: `col_0based = h3_pos + 16 - 1` (e.g. F159 → col 174).
- Proof: `COORD_VALIDATION.json` / `scripts/validate_ha_coords.py` (TX/50 WT `data/h3n2/wt_ha.txt` + RBS landmarks Y98/W153/H183/Y195).
- Anc AA under `data/train` are **ungapped** (`total_gap_chars=0`); pad-to-566 ≠ homology.
- **H1** uses different numbering + ~17-aa signal — do **not** reuse H3 offsets; H1 sites listed for transfer only.

## Curated H3 site sets (flu only)

**Primary (train lit mask; unique H3 positions):**  
145, 155, 156, 158, 159, 183, 189, 193, 225  
(= H155T/Q156H SY97→FU02; F159S/Y 2014–15; H183L; K189E; N225D; Fig. 5 Koel-colored 145/158/193).

**Secondary (Fig. 6 TX/50 escape genotypes):**  
75, 88, 94, 107, 122, 127, 128, 144, 157, 160, 172, 174, 192, 197, 203, 207, 217, 219, 220, 242, 246.

**Region prior:** globular-head library span H3 **#63–252** → `mut_hotspot_mask_h3_globular_head.pt`.

**H1 transfer (not in H3 mask):** Sa **153–156** (+ D127E); broader Sa/Koel set in `flu_lit_sites.py`. Needs separate H1 coord validation before train.

## Experimental / modeling claims (from author MS)

- Escape screens: random globular-head libraries × human/ferret convalescent sera; HI + antigenic cartography.
- H1: Sa **153–156** dominate antigenic shift; single AA can give ≥16-fold HI drop; field variants match lab escapes; ferret/mouse challenge shows immune evasion.
- H3 retrospective: SY97→FU02 needs **H155T+Q156H**; intermediate genetic background matters (Kwangju/219 already 155T).
- H3 contemporary (TX/50): HA-**159** most frequent escape site; F159Y/S map toward 3C.2a/3C.3a; **K189E+N225D** can advance antigenically without 159Y; method predicts **cluster location**, not exact next sequence; NA can matter; **re-run as clades update**.
- Assays: HI (± oseltamivir for H3), cartography, mouse/ferret challenge — not DMS fitness tables.

## Local MSA vs lit (`data/train`, L=566)

Against 9 unique primary H3 positions, MSA mut-freq top-15% (85/566) has **high recall, low precision**:

| selector | precision | recall | missed |
|---|---|---|---|
| top-9 | 0.22 | 0.22 | most |
| top-32 | 0.19 | 0.67 | … |
| top-85 (≈frac 0.15) | ~0.10 | **0.89** | **H183** only (rank ~155) |

Prefer **`mut_hotspot_mask_nmicrobiol_lit.pt`** for antigenic mut recovery (keeps H183); MSA mut-freq remains diagnostic / soft prior. (Unlike COVID geo MSA, flu polymorphism already covers most lit sites — still use lit for narrative KPI + force.)

## Actionable train / eval (H3N2 only)

1. **Hotspot:** `--mut-hotspot-mask results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt` + `--mut-hotspot-weight 5` + prefer `--mut-hotspot-force` (escape sites can look “cons” on a bridge step).
2. **Optional OR:** `mut_hotspot_mask_flu_lit_or_head.pt` if recall of rare head sites matters.
3. **Eval KPI:** `--lit-hotspot-mask results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt` → `flu_hotspot_mut_frac` / `gt_flu_hotspot_mut_frac`.
4. **Do not** point flu jobs at COVID PMC/RBD masks (refused by eval guard).
5. **Temporal:** paper is about *future* antigenic clusters → prefer `data/h3n2` temporal / forecast splits over mixed `data/train`.
6. **Fitness vs escape:** N225D is RBS/avidity as well as antigenic; K189E antigenic. Dual masks or higher weight on {155,156,159,189} vs RBS-only if narrating Table-style escape.
7. **Epistasis:** SY97→FU02 and Fig. 6 multi-mut escapes → eval joint recovery on {155,156}, {159,225}, {189,225}, not only marginals.
8. **Background dependence:** intermediate genotypes matter → temporal / clade-holdout > i.i.d. leaf holdout for mut_recovery claims.
9. **Reference process:** soft empirical entropy + hard lit mask mirrors “serum-selected escape sites”; keep `--entropy-source empirical`.
10. **H1:** separate mask + signal validation; do not transfer H3 columns to `data/h1n1`.

## Modeling brainstorm (concrete flags / experiments)

Stack assumptions: TreeSBM bridge, `log R_θ = log R0 + c_θ`, optional semigroup `λ_semi`, fitness-β R0 tilt, R0 backends, hotspot force/weight, mrs sweeps, deep mut head, temporal splits.

| # | Experiment | Flags / paths | Why (paper ↔ stack) |
|---|---|---|---|
| A | **h3n2_v*_lit_hotspot** (primary) | v5-style mutrec + `MUT_HOTSPOT_MASK=…/mut_hotspot_mask_nmicrobiol_lit.pt` `WEIGHT=5` `FORCE=1` `λ_semi=0.05` | Curated antigenic sites → L_mut; force keeps 183/159 learnable when `aa_t==x1` |
| B | lit ∪ head band | `mut_hotspot_mask_flu_lit_or_head.pt` FORCE=1 | Library span 63–252 = soft region prior like COVID RBD band |
| C | MSA mut_freq only | `MUT_HOTSPOT_SCORE=mut_freq` `FRAC=0.15` (no lit mask) | Ablation; expect worse H183 / narrative KPI |
| D | Fitness-β × lit | A + `FITNESS_BETA∈{0.5,1.0}` `FITNESS_SCORE=log_R0` | Escape ≠ growth; tilt R0 while lit hard-mask handles antigenic columns |
| E | R0 backend sweep | `--r0-backend esm2` vs larger/stub (match precompute tag) under A | Does better R0 reduce need for FORCE on conserved antigenic sites? |
| F | Temporal lit | Train `data/h3n2/train` (or forecast cut) + eval future leaves; report `flu_hotspot_mut_frac` @ mrs 0.3/0.5 | Paper = predict next cluster before surveillance dominance |
| G | mrs × site-softmax | Eval A @ mrs∈{0.3,0.5,1.0} ± `--site-softmax-sample` | Mut recovery ops; lit frac should rise with mrs if hotspot learning worked |
| H | Deep mut head | A + `DEEP_MUT=1` | More capacity for multi-site Sa/B patterns (159+189+225) |
| I | Bundle eval | Post-hoc joint hit-rate on {155,156}, {189,225}, {159,Y/S} | Paper epistasis / cluster motifs |
| J | H1 transfer (later) | New H1 mask + `--signal-len` validation; **never** H3 `.pt` on H1 | Sa 153–156 is H1-native story |

**Architecture interactions:**

- Hotspot **weight** stacks with entropy soft weights; **force** removes hotspot cols from `cons_mask` — use with small lit sets (n=9) so cons retention elsewhere stays healthy.
- Semigroup `L_semi` regularizes rate composition across times; orthogonal to site mask but helps when FORCE expands mut sites.
- Fitness-β tilts **R0** before `c_θ`; lit mask tilts **loss** — complementary, not substitutes.
- PSSM gate / mut-aa-emb (v5) help AA identity at known escape sites once site is selected.
- Do not bake HI cartography into R0; keep lit mask as external evaluator + train prior (same philosophy as COVID EVEscape-as-eval).

## CLI

```bash
# Validate H3 coords + (re)write flu lit / head masks:
python scripts/validate_ha_coords.py \
  --data data/train --max-seq-len 566 \
  --out-dir results/flu_mutfreq_vs_lit

# MSA mut-freq vs lit (prefer data/h3n2/train on Betty when trees exist):
python scripts/compute_flu_mut_freq_vs_lit.py \
  --data data/train --max-seq-len 566 \
  --out-dir results/flu_mutfreq_vs_lit --hotspot-frac 0.15

# Enrichment KPI (explicit flu mask only):
python scripts/eval_evescape_enrichment.py \
  --checkpoint checkpoints/h3n2_v2/best.pt \
  --data data/h3n2/test --max-seq-len 566 \
  --lit-hotspot-mask results/flu_mutfreq_vs_lit/mut_hotspot_mask_nmicrobiol_lit.pt \
  --mutation-rate-scale 0.3 \
  --out checkpoints/eval_enrichment_h3n2_flu_hotspot_mrs0.3.json
```

Betty train paste (prepared, **not submitted**): `scripts/_paste_betty_h3n2_lit_hotspot.sh`.
