# ICLR submission: exact numbers to update

Scope: `TreeSBM__ICLR_ (10).pdf`. Companion to `SPIKE_FRAME_BUG.md` (cause) and
`NUMBERS_THAT_CHANGE.md` (recomputed VOC-panel values).

The headline: **every number in the paper that is computed by comparing
sequences to each other survives. Every number computed by comparing residues at
a shared column index has to be recomputed**, and for SARS-CoV-2 and HIV the
recomputation is not cosmetic.

The good news is that it is cheap. Generation caches already exist on Betty for
essentially every table, so this is a CPU rescoring pass, not a GPU regeneration
pass. See §5.

---

## 1. Why the split falls where it does

`benchmarks/metrics/sequences.py::positional_recovery` is the function behind
mutation recall, site recall, `AA|hit`, and conservation retention:

```242:242:benchmarks/metrics/sequences.py
    L = min(len(root), len(gt), len(gen))
```

It then reads `root[i], gt[i], gen[i]` at the same `i`. That is only meaningful
if the root, the observed leaf, and the generated leaf are in one frame.

Two facts about the pipeline decide when they are:

1. **TreeSBM only substitutes** ("TreeSBM models single substitution mutations
   only", §4.4). A generated leaf therefore always inherits the root's exact
   length and frame. Generated-vs-root is always aligned.
2. **Observed leaves carry their own indels.** Observed-vs-root is aligned only
   when the leaf happens to have the same deletion pattern as the inferred
   ancestral root.

So the failure mode is specifically **observed-leaf-vs-root**, and separately,
any mask defined in *reference* coordinates is wrong whenever the **root** is
off-reference, because the whole generated tree sits in the root's frame.

### Exposure per dataset

`scripts/audit_frame_all_datasets.py`, output in `frame_audit_all_datasets.json`.

| Dataset | trees | modal len | roots at modal | leaves at modal | leaves in root's frame |
|---|---|---|---|---|---|
| SARS-CoV-2 Spike (Brazil test) | 40 | 1271 | 0.225 | 0.465 | **0.467** |
| SARS-CoV-2 Spike (train) | 335 | 1271 | 0.164 | 0.326 | **0.477** |
| Influenza H3N2 HA | 45 | 566 | 0.889 | 0.914 | **0.936** |
| Influenza H1N1 HA | 40 | 566 | 0.725 | 0.872 | **0.744** |
| HIV Env (geo test) | 14 | 853 | 0.071 | 0.087 | **0.018** |

Read this as the fraction of leaf comparisons that are valid today:

- **HIV Env: 1.8%.** Env's variable loops differ in length between essentially
  every isolate, so almost no leaf shares its root's frame. Every column-wise
  HIV number in the paper is noise.
- **SARS-CoV-2: 47%.** Roughly half the comparisons are misaligned, and only
  22.5% of test roots are even at the modal length — so reference-coordinate
  masks land on the wrong residues for 77.5% of roots.
- **H1N1: 74%.** Real but moderate.
- **H3N2: 94%.** Small; expect third-decimal movement.

---

## 2. Cell-by-cell change list

Legend — **SAFE**: no reference indexing, do not touch. **RESCORE**: column-wise
or mask-indexed, must be recomputed. **CLAIM**: wrong independent of the bug.

### Table 2 — H3N2 tree topology (p. 9)

RF, quartet, branch W1, terminal edit distance are all sequence-to-sequence or
topology comparisons. **SAFE, entire table.**

### Table 3 — H3N2 future lineage forecasting (p. 10)

| Column | Values | Status |
|---|---|---|
| Coverage@10 / @100 | 0.662–0.725 | SAFE |
| Median min. dist | 2.65 / 2.65 / 2.40 / 2.24 | SAFE |
| **Mutation recall** | 0.020 / 0.004 / 0.079 / **0.156** | RESCORE (H3N2, ~6% of comparisons affected — expect small moves) |
| Clade recall | 0.000 / 0.000 / 0.300 / 0.100 | RESCORE if clades are defined by reference-position mutation sets |

Also: TreeSBM's clade recall (0.100) is **worse than the autoregressive
baseline's (0.300)** in the paper's own headline forecasting table, and the text
does not acknowledge it. Decide whether to drop the column or address it.

### Table 4 — Viral evolution forecasting (p. 11). Highest priority.

**SARS-CoV-2 Spike** — 53% of leaf comparisons misaligned, 77.5% of roots
off-reference.

| Metric | PLM / AR / TreeSBM | Status |
|---|---|---|
| Coverage@10, Coverage@100 | 0.775 / 0.775 / 0.788 | SAFE |
| Min. distance | 1.425 / 1.425 / 1.388 | SAFE |
| **Key mutation recall** | 0.048 / 0.006 / **0.219** | RESCORE — this is the paper's central quantitative claim |
| **AA \| hit** | 0.017 / 0.094 / **0.375** | RESCORE |
| **Conservation retention** | 0.985 / 0.993 / **0.840** | RESCORE |
| **Antigenic recall** | 0.00 / 0.017 / **0.005** | RESCORE (mask in reference coords) |
| **EVEScape Δ** | -0.424 / -0.371 / **-0.048** | RESCORE (score tensor in reference coords) |

**Influenza H1N1 HA** — 26% misaligned, 27.5% of roots off-reference. Same five
metrics RESCORE: key mutation recall (0.071 / 0.006 / 0.241), AA|hit (0.136 /
0.008 / 0.618), conservation retention, antigenic recall (0.115 / 0.00 / 0.163),
EVEScape Δ. Coverage and min. distance SAFE.

Two CLAIM problems in this table, unrelated to the bug:

- **Coverage@10 equals Coverage@100 in every SARS-CoV-2 row.** Ten-fold more
  samples buying exactly zero additional coverage says the metric is saturated
  or K is not actually varying. A reviewer will read this as a bug. Diagnose
  before resubmission.
- §5.4 says TreeSBM keeps "negligible drops in fitness," but its conservation
  retention (0.840) is the **lowest in the table** by ~15 points. Either the
  sentence goes or the metric does.

### Table 5 — Reference process ablation (p. 12)

| Column | Status |
|---|---|
| Tree-KL (0.693 × 7), Coverage@100 (0.775 × 7) | SAFE — but see below |
| pLM NLL | SAFE |
| **Mut. recall** (0.023 / 0.018 / 0.001 / 0.035 / 0.036 / 0.012 / **0.038**) | RESCORE |

**CLAIM:** Tree-KL is 0.693 and Coverage@100 is 0.775 for *all seven rows*, and
the same two constants reappear in Tables 6, C.1, C.2, D.1, and E.2. 0.693 is
ln 2, i.e. the Tree-KL estimator is returning its uninformative value everywhere.
This is the single most damaging thing in the submission: it says two of the
three headline metrics are insensitive to every ablation the paper runs. Fix the
estimator or remove both columns — rescoring will not help.

### Table 6 — Core TreeSBM ablations (p. 13)

RESCORE the **Cons. retention**, **AA | hit**, and **Mut. recall** columns for all
six rows. Tree-KL, Coverage@100, pLM NLL as per Table 5.

Note `AA|hit` is identical (0.375) for "without branch-length head" and "full",
as are Mut. recall (0.08) and pLM NLL (0.457) — the two rows are the same run.
`coverage_cache_covid_t8_no_bl_eval` exists separately, so rescore both properly.

Also: §5.6 says the table removes "internal-node sequence modeling," but no such
row exists. The cache `coverage_cache_covid_t8_no_internal_seqs` is on Betty, so
add the row.

### Table 7 — Evolutionary panels (p. 13)

Entirely empty. Not a correction; see §4.

### Table C.1 — Full tree generation by dataset (p. 19)

Tree-KL, Split-KL, RF, quartet, branch W1, pLM NLL. **SAFE, entire table**
(modulo the Tree-KL estimator problem above).

### Table C.2 — Full future-lineage forecasting (p. 20)

Coverage and min. distance SAFE. RESCORE the **Mut. recall** and **Clade recall**
columns for all four datasets. The HIV Env values (0.015 / 0.026) are computed on
1.8% valid comparisons and should be treated as unreported until rescored.

**CLAIM:** Table 3 gives H3N2 TreeSBM mutation recall as **0.156**; Table C.2
gives the same quantity as **0.100**. Coverage and min. distance agree across the
two tables, so this is a genuine inconsistency, not a definition difference.
Likewise Table 4's SARS-CoV-2 "key mutation recall" (0.219) versus Table C.2's
"Mut. recall" (0.080) — if these are different metrics, define both.

### Table C.3 — Forecast horizon (p. 20)

RESCORE the Mut. recall column (H3N2 0.099 / 0.131; HIV 0.001 / 0.050). Coverage
and min. distance SAFE. Most cells are dashes; caches `coverage_cache_c3_*`
exist for covid, h1n1, h3n2, hiv_geo, hiv_temporal at N16, so the table can be
completed in the same pass.

### Table D.1 — Reference process components (p. 20)

RESCORE the **AA | hit** column (0.108 / 0.098 / 0.375 / 0.086 / 0.375).
Coverage@100 and Tree-KL as per Table 5.

Note "pLM + branching" and "TreeSBM full bridge" both report 0.375, and
"pLM + fitness + branching" reports 0.086 — i.e. adding fitness on top of
branching *destroys* AA|hit while the full bridge recovers it exactly. That
pattern needs an explanation or it reads as a copy error.

### Figures

| Figure | Status |
|---|---|
| Fig. 3 — coverage vs ε | **SAFE.** Hamming-threshold curve, no reference indexing. |
| Fig. 4 — Gamma trajectory (MZ397166.1) | **SAFE, verified.** P.1 carries no Spike deletion; g003's root is 1273 aa and 98.7% of its leaves share that frame. Recomputed `bundle_frac` 0.7133 → 0.7167. |
| Fig. 5 — Gamma escape recovery, n=300 | **SAFE.** Generated leaves inherit the 1273-aa root frame, so the reference-position mask and EVEscape lookup are correctly placed. The "eleven of twelve / four in more than half" claim stands. |
| Fig. 2, H.1, H.2 | SAFE (visualizations). |

Gamma is the control for the whole audit: it has no Spike deletion, so if the
alignment fix were wrong, Gamma would move. It does not.

### Table J.1 — Linear probes (p. 24)

Not affected by the frame bug, but **CLAIM**: the trained graph transformer
(GraphTF) is the *worst* row for AncAA (0.086 vs raw ESM's 0.748) and MutPos
(0.009 vs 0.396), and has RootDist R² of -13.16. The table is captioned as
validation but shows the encoder discarding the sequence information the
mutation head depends on. Either reframe it as a diagnostic or cut it.

### §5.3 text (p. 9)

> "With 100 generated trees, 72.5 percent of generated leaves were within two
> mutations of observed 2024 sequences"

Backwards. Table 3's own caption defines Coverage@K as the proportion of
*observed sequences* for which at least one *generated leaf* lies within ε — not
the fraction of generated leaves. Rewrite; this one is free and a reviewer will
catch it.

---

## 3. Also worth stating accurately: training-set VOC prevalence

From `train_leakage_audit_fixed.json` (recomputed with the alignment fix), the
defining bundles are far more common in training than the pre-fix audit showed:
Delta's signature appears in **47% of all training leaves**, Omicron BA.1's in
38%, Alpha's in 36%.

This does not break the evaluation — the standard we apply is an unseen *root*
and unseen *tree*, and shared mutations across a pandemic corpus are expected.
But it does mean the paper cannot imply these signatures were rare in training.
State the tree-level standard explicitly in §5.4 or the appendix.

---

## 4. Empty tables that the existing caches can already fill

These are dashes in the submission, and Betty already has the generations:

| Table | Cache on Betty |
|---|---|
| E.1 Bridge matching ablations | `coverage_cache_covid_e1_e1_no_doob`, `_no_terminal`, `_terminal_only` |
| E.2 Architecture ablations | `coverage_cache_covid_e2_e2_no_mut_head`, `_no_stop_head` (+ t8 caches) |
| E.3 Sampling sensitivity | `coverage_cache_covid_e3_temp_0p5`, `_temp_1p5`, `_maxl_100`, `_maxl_800`; K=500 caches exist for covid/h1n1/h3n2/hiv |
| C.3 Horizon | `coverage_cache_c3_{covid,h1n1,h3n2,hiv_geo,hiv_temporal}_N16` |

Tables 7 (panels), F.1 (data scaling), F.2 (cross-domain), G.1 (likelihood
ranking), I.1 (cost), J.2/J.3 (hyperparameters) have no caches. J.2 and J.3 are
just transcription from configs and should be filled regardless — an ICLR
submission with an empty hyperparameter table is an easy desk reject.

---

## 5. Execution plan

The fix lives in the metric, not the model, and TreeSBM's generations are cached,
so **nothing needs to be regenerated on GPU.**

`benchmarks/coverage_curves.py` already supports `--rescore-from DIR`, which
loads `tree_*.json` and skips generation. The work is therefore:

**Step 1 — make the metric alignment-aware (½ day).**
Give `positional_recovery` the same treatment `voc_threat_lib.seq_at` got: map
root, observed leaf, and generated leaf onto reference coordinates via
`ref_pos_map` before comparing, and score only positions present in all three.
Sequences already at reference length skip alignment, so H3N2 stays fast.
Apply the same mapping to the mask lookups in `score_cache_antigenic.py` and
`eval_evescape_enrichment.py`.

**Step 2 — regression-test on Gamma (½ day).**
Re-run Fig. 4 / Fig. 5 numbers. They must not move. If they do, the fix is wrong.

**Step 3 — rescore, CPU only (1 day wall-clock).**
Sweep `--rescore-from` over the caches in the order the tables matter:
Table 4 (`*_viral_eval_K100`) → Tables 3 and C.2 (`*_N16`) → Table 5 (`t7_*`) →
Tables 6 and D.1 (`t8_*`) → Tables C.3 (`c3_*`), E.1–E.3 (`e1_*`, `e2_*`, `e3_*`).

**Step 4 — the two things rescoring cannot fix (1–2 days).**
- Tree-KL pinned at ln 2 and Coverage@100 pinned at 0.775 across every ablation.
- Coverage@10 == Coverage@100 on SARS-CoV-2.
Both are estimator or harness bugs and both are more dangerous to the paper than
the frame bug, because they make the ablations vacuous rather than merely wrong.

**Step 5 — retraining, only if needed (open).**
The `covid_v7_pmc_hotspot` checkpoint uses a reference-coordinate literature mask
*during training*, applied to padded columns. For the 77.5% of roots that are
off-reference the mask supervises the wrong residues. Steps 1–3 correct the
*reporting*; correcting the *training* means a retrain and a re-run of the
`t8_no_lit_mask` comparison. Decide after Step 3 shows how large the reporting
correction is.

---

## 6. Added after Step 1: truncated ancestral roots

Found while regression-testing the fix. **28% of COVID trees (train and test) and
64% of HIV trees have a root that is a short N-terminal fragment** — some as short
as 6–27 aa — while their leaves are full length. Because generation starts at the
root and only substitutes, every generated leaf in those trees is equally
truncated, so they structurally cannot contain the RBD.

This is worse than the frame bug in two ways. It affects **training**, not only
reporting (94 of 335 COVID training trees). And because `hamming` and `identity`
truncate to the shorter sequence, it also reaches **Coverage@K and min-distance**
— the metrics marked SAFE above are safe from the frame bug but not from this.

Details and proposed fix in [`TRUNCATED_ROOTS.md`](TRUNCATED_ROOTS.md). It has to
be settled before the rescoring pass, or the corrected numbers get recomputed on
the same broken subset.

## 7. Expected direction of the correction

The regression test measured legacy vs aligned scoring on an off-frame tree
(g040, 2% of leaves in frame):

| metric | legacy → aligned |
|---|---|
| `mut_recovery` | 0.0760 → 0.0157 |
| `cons_retention` | 0.8468 → 0.7111 |
| `aa_acc_given_hit` | 0.2524 → 0.0411 |

A frame shift makes nearly every position look like a mutating site, so the model
collects chance-level credit across thousands of spurious ones. Correcting it
shrinks the denominator to the real mutations, where performance is much lower.

**Expect Table 4's SARS-CoV-2 `AA|hit` (0.375) and key mutation recall (0.219) to
fall.** Plan the narrative accordingly rather than being surprised by it.

## 8. Priority order

1. **Truncated roots** (§6) — rebuild or exclude. Everything downstream is
   computed on this corpus, so it goes first.
2. Step 4's Tree-KL / coverage-saturation diagnosis — vacuous ablations sink the
   paper faster than wrong ones.
3. Table 4 SARS-CoV-2 rescore — the central claim.
4. §5.3 coverage sentence and §5.4 "negligible fitness drop" sentence — free.
5. Tables 3, 5, 6, C.2, D.1 rescore.
6. HIV: either rescore and report honestly, or withdraw the column-wise HIV
   numbers. At 1.8% valid comparisons the current values are indefensible.
7. Fill J.2/J.3, then E.1–E.3 and C.3 from cache.
8. Table J.1 reframing.
