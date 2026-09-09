# Appendix E.3 — Sampling hyperparameter sensitivity

**Status 2026-08-17:** COVID one-factor grid **temp / max_leaves COMPLETED** (**7639276–79**). K=50/100 from existing Table5 eabs; K=500 pending C.2 **7635911** RUNNING (CSV 0 bytes).

**Paper columns:** Setting, Value, Coverage, Diversity, Tree-KL, Runtime  
**Additive (on top):** Cons, AA|hit, Mut (when in CSV). Antigenic **—** (coverage runner lacks lit hotspot frac).

**Protocol:** COVID Brazil geo-test; ckpt `covid_v5_mutrec`; same 5 Table5 roots; N=16; primary = `coverage_obs_e2`.  
**Design:** one-factor-at-a-time (not 3³ factorial). Defaults held: K_max=100, `site_temperature=1.0`, max_leaves = production N-adapter (documented as 400 gen default).

**Grid:**
| Axis | Values | Jobs / source |
|---|---|---|
| K | 50, 100, 500 | existing T5 eabs; K500 ← **7635911** |
| temperature | 0.5, 1.0, 1.5 | **7639276** / baseline / **7639277** |
| max_leaves | 100, 400, 800 | **7639278** / production default / **7639279** |

Job IDs: [`table_e3_job_ids.json`](table_e3_job_ids.json). Submit: `scripts/betty_submit_e3.sh`.

**Sources (mapped):**
- H3N2 K-sweep ε=5: `table4_coverage_vs_K_eps5.csv` / `coverage_curves_h3n2_N16_eabs.csv` (job **7459412**)
- COVID K-sweep e-abs: `coverage_curves_covid_N16_table5_eabs.csv`
- Runtime = `runtime_gen_sec` for the **K_max pool** (same wall for all K rows — pool generated once, then scored at each K)

---

## Paste-ready

| Setting | Value | Cov ↑ | Diversity ↑ | Tree-KL ↓ | Runtime (s) | Cons ↑ | AA\|hit ↑ | Mut ↑ | Status |
|---|---|---:|---|---|---:|---:|---:|---:|---|
| K (H3N2, Cov@e5) | 50 | **0.812** | — | — | 15668† | **1.000** | — | **0.252** | DONE (existing curve) |
| K (H3N2, Cov@e5) | 100 | **0.812** | — | — | 15668† | **1.000** | — | **0.156** | DONE |
| K (H3N2, Cov@e2) | 50 | **0.738** | — | — | 15668† | **1.000** | — | **0.252** | DONE |
| K (H3N2, Cov@e2) | 100 | **0.725** | — | — | 15668† | **1.000** | — | **0.156** | DONE |
| K (COVID, Cov@e2) | 50 | **0.775** | — | — | 21421† | **1.000**‡ | **1.000**‡ | 0.001 | DONE |
| K (COVID, Cov@e2) | 100 | **0.788** | — | — | 21421† | **1.000**‡ | **1.000**‡ | 0.011 | DONE |
| K (COVID, Cov@e2) | 500 | — | — | — | — | — | — | — | **RUNNING** C.2 **7635911** (CSV 0 bytes) |
| Temperature | 0.5 | **0.775** | — | — | 20498 | **1.000**‡ | **1.000**‡ | 0.001 | **7639276** COMPLETED |
| Temperature | 1.0 | **0.788** | — | — | 21421† | **1.000**‡ | — | 0.011 | DONE (= Table5 treesbm @K=100) |
| Temperature | 1.5 | **0.775** | — | — | 20639 | **1.000**‡ | **1.000**‡ | 0.001 | **7639277** COMPLETED |
| Max leaves | 100 | **0.775** | — | — | 24655 | **1.000**‡ | **1.000**‡ | 0.039 | **7639278** COMPLETED |
| Max leaves | 400 | **0.788** | — | — | 21421† | — | — | — | DONE (production default / Table5) |
| Max leaves | 800 | **0.775** | — | — | 44789 | **1.000**‡ | **1.000**‡ | 0.042 | **7639279** COMPLETED |

† Wall seconds for generating the K≤100 pool over **5 roots** (not per-tree). COVID ≈ **5.95 h**; H3N2 ≈ **4.35 h**.  
‡ Coverage-path cons/aa|hit on near-identical Brazil leaves — **not** comparable to enrich mrs0.5 additives in D.1/T7. Prefer enrich for sequence quality.

**Diversity:** paper column stays **—**. Coverage CSV has no dedicated diversity (e.g. mean pairwise leaf Hamming among gens). Do not invent from `frac_gen_*`.

**Tree-KL:** **—** on all E.3 rows (no per-K Tree-KL pass; optional later).

---

## Caption honesty

One-factor COVID TreeSBM sensitivity around production defaults, plus mapped H3/COVID K-curves. New T≠1 / max_leaves≠400 jobs (**7639276–79**) scored Cov@e2=**0.775** on the same 5 Brazil roots (flatness; Table5 T=1.0 / maxl≈400 mapped row stays **0.788**). Rows without artifacts stay **—**. Full factorial T×K×max_leaves **not** run.
