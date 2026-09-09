# PMC10142771 → TreeSBM modeling notes

Source: Kumar et al., *Viruses* 2023, PMC10142771.  
Coords: **1-based UniProt P0DTC2 / Wuhan-Hu-1 spike**. Our MSA columns are **0-based**; ungapped Wuhan-frame AA ⇒ `spike_pos = col + 1`.

Paper site selection is **not** “any polymorphic column”: mutations were chosen by **global entropic score + emergence/spread/transmission + RBD location**, with **D614G as universal reference**. That matches our MSA-select → tree-apply design, but raw MSA mut-freq alone can miss rare antigenic RBD sites — use curated lit / RBD masks when recall vs primary set is low.

## Curated site sets (use for hotspot validation)

**Primary (study-selected + D614G + S477N):**  
A222V, L212I (NTD); Y369C, K417N, N439K, L452R, Y453F, S477N, T478K, E484K, G496S, N501Y (RBD/RBM); T547K, D614G (S1); N856K, L981F (S2).

**Secondary (VoC body text):** Omicron RBD cluster G339D, S371L/S373P/S375F, N440K, G446S, E484A, Q493R, Q498R, Y505H; P681R; K417T / E484Q.

**Regions:** RBD ~319–541 (paper cites NTD–RBD contact near 333–528); RBM loop **S477–G485** called highly flexible / ACE2-contact; interface residues Y369, N370, S383, K417, Y421 flagged mutation-prone by SASA.

Artifacts: `mut_hotspot_mask_pmc_lit.pt`, `mut_hotspot_mask_rbd_band.pt`, `precision_recall_vs_lit.csv`.

## Actionable train / eval changes

1. **Hotspot score:** prefer `--mut-hotspot-score mut_freq` (frac ≠ consensus) or load **`mut_hotspot_mask_pmc_lit.pt`** when MSA top-15% recall vs primary lit &lt; ~0.75. Keep entropy as soft weighting; hard mask = antigenic sites.
2. **Region prior:** optional OR of lit mask with RBD band (319–541) so RBM escape sites stay in `L_mut` even if locally fixed in train MSA.
3. **Eval KPI split:** report `mut_recovery` / `site_recall` / `aa_acc_given_hit` **restricted to primary lit columns** and to RBD band — paper cares about antigenic/RBD recovery, not uniform column average.
4. **EVEscape / enrichment:** already RBD-focused; align enrichment site set to primary+Omicron RBD list above (esp. 417/484/501/452/478/477).
5. **Fitness vs escape:** D614G = fitness/transmission; E484K/N439K/Y453F/L452R = escape. Do **not** collapse into one hotspot weight — consider dual masks (fitness vs escape) or higher weight on escape RBD sites for Table 5 narrative.
6. **Epistasis / co-occurrence:** paper repeatedly notes linked sets (N501Y+N439K+Δ69/70; Beta 417+484+501; Delta 452+478; Omicron multi-RBD). Tree bridge is edge-local; add **eval** of joint recovery on known VoC mutation bundles, not only marginal site metrics.
7. **Temporal / geographic:** Nextstrain-style waves matter — prefer `data/covid_temporal` / clade-holdout for mut_recovery so train MSA hotspots match future leaves (geo MSA can overfit shared D614G-era polymorphisms).
8. **Reference process:** soft floor+αH on empirical entropy is closest to their “global entropic score”; hard `mut_freq` or lit mask is the discrete analogue of “selected mutations.” Keep `--entropy-source empirical` for soft weights even when hard score is `mut_freq`.
9. **Cons vs mut:** RBD “interface / open-state” sites can look conserved on a given bridge step (`aa_t == x1`) yet be antigenically critical — `--mut-hotspot-force` on lit/RBD mask is justified for those columns.
10. **S2 / stability sites** (L981F, N856K, T547K): include in lit mask but expect lower EVEscape signal; don’t use EVEscape alone to define hotspots.

## Validation (COVID train on Betty: 335 trees, 146495 node seqs)

Against PMC primary **16 unique** spike positions, MSA mut-freq top-15% (192/1280) has **poor lit recall**:

| selector | precision | recall | hits |
|---|---|---|---|
| top-16 MSA mut-freq | 0.00 | 0.00 | — |
| top-64 | 0.016 | 0.062 | Y369C |
| top-96 / top-192 (≈frac 0.15) | ~0.01–0.02 | **0.125** | L212I, Y369C only |

Many VoC RBD sites are mid-ranked by raw polymorphism on this geo train MSA, so **prefer `mut_hotspot_mask_pmc_lit.pt` (or lit ∪ RBD band)** for antigenic mut recovery; keep MSA mut-freq as diagnostic / soft prior.

## Coordinate validation (Wuhan / P0DTC2)

### Q: Normalize length first to find the spike window?

**No — same length alone is not enough.** Pad/truncate to `max_seq_len=1280` does **not** make column `j` homologous across sequences. You need a **column-aligned MSA** (or map each sequence to P0DTC2). After that, shared width `L` is a spike window; it is 1:1 with P0DTC2 only if the alignment is **reference-anchored**. Landmark AA checks (D614, N501, …) are the proof.

Empirical check of `group_*_anc_aa.fasta`: **zero `'-'` gap characters** → left-aligned ungapped AA strings, not a gapped MSA. Many trees mix lengths; absolute columns are Wuhan-frame only for full-length `L=1273` (and same indel pattern).

**Mapping is correct** for full-length ungapped spike (`L=1273`): `spike_pos = col + 1`.
Proof: `COORD_VALIDATION.json` (`scripts/validate_spike_coords.py`).

| check | result |
|---|---|
| `data/covid/wt.txt` (P0DTC2) | D@614, N@501, K@417, E@484, … |
| EVEscape `reference_seq` | identity 1.0 to WT; scored RBD 331–531 |
| Train **L=1273 only** (~38k seqs) | mean id ≈ 0.99; cons@614=**G** (D614G), cons@501=**N**, cons@417=**K** |
| Motif `YQGVNCT` | starts at col **611** on L=1273 (= residues 612–618) |
| Gap / MSA audit | `total_gap_chars=0`; pad≠homology |

Caveat: shorter seqs (Δ69/70 etc.) without gaps shift absolute columns (e.g. L=1271 `YQGVNCT` at 609, not 611). PMC/EVEscape masks stay Wuhan-frame (exact on L=1273). Earlier `msa_mut_freq_by_site.csv` consensus mixed lengths — ignore for indexing.

Eval KPI: `pmc_hotspot_mut_frac` / `gt_pmc_hotspot_mut_frac` in `eval_evescape_enrichment.py` JSON
(= mean over leaves of \|mut_cols ∩ PMC\| / \|mut_cols\|; gen≠root vs GT≠root).

### Leaf-only vs tree-wide metrics

**Primary** `mut_recovery` / `site_recall` / `aa_acc_given_hit` / `cons_retention` are **leaf-only**:
best-match generated **leaf** vs GT **leaf** vs root. Internal/ancestral nodes are not scored.
This is intentional for tip forecasting. Timing of when a mutation appears on an internal edge
does not enter the primary KPI — only the leaf AA does.

Enrichment JSON also emits **extra** (non-replacing) companions:
- `mut_recovery_any_descendant` / `site_recall_any_descendant` — credit if **any** gen leaf recovers a GT leaf mut
- `mut_site_recall_path_union` — union of mut sites over GT leaves vs gen leaves

### Domain mut fractions (COVID spike)

Alongside `pmc_hotspot_mut_frac`, enrichment (auto when `max_seq_len≥1000`) emits:

| key | band (1-based Wuhan/P0DTC2) |
|-----|------------------------------|
| `ntd_mut_frac` / `gt_*` | NTD 13–305 |
| `rbd_mut_frac` / `gt_*` | RBD 319–541 |
| `rbm_mut_frac` / `gt_*` | RBM 438–506 |
| `s1_mut_frac` / `gt_*` | S1 13–685 |
| `s2_mut_frac` / `gt_*` | S2 686–1273 |

Same definition as PMC hotspot frac (leaf mut columns). Refused for flu lengths. See `src/spike_domains.py`.

## CLI

```bash
# Validate coords + regenerate PMC / RBD masks (Betty):
python scripts/validate_spike_coords.py \
  --data data/covid/train --max-seq-len 1280 \
  --out-dir results/covid_mutfreq_vs_lit

# Recompute MSA mut-freq vs lit (Betty; trees present):
python scripts/compute_msa_mut_freq_vs_lit.py \
  --data data/covid/train --max-seq-len 1280 \
  --out-dir results/covid_mutfreq_vs_lit --hotspot-frac 0.15

# Preferred antigenic train (curated PMC mask → tree-apply):
HOTSPOT=1 CKPT_DIR=checkpoints/covid_v7_pmc_hotspot \
  MUT_HOTSPOT_MASK=results/covid_mutfreq_vs_lit/mut_hotspot_mask_pmc_lit.pt \
  MUT_HOTSPOT_FRAC= MUT_HOTSPOT_TOPK= \
  MUT_HOTSPOT_WEIGHT=5 MUT_HOTSPOT_FORCE=1 \
  sbatch --qos=mig-max --export=ALL,HOTSPOT,CKPT_DIR,MUT_HOTSPOT_MASK,MUT_HOTSPOT_FRAC,MUT_HOTSPOT_TOPK,MUT_HOTSPOT_WEIGHT,MUT_HOTSPOT_FORCE \
    scripts/slurm_covid_train_mut_recovery.sh

# Or MSA mut-freq hard hotspots:
HOTSPOT=1 MUT_HOTSPOT_SCORE=mut_freq MUT_HOTSPOT_FRAC=0.15 \
  CKPT_DIR=checkpoints/covid_v6_mutfreq \
  sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh
```
