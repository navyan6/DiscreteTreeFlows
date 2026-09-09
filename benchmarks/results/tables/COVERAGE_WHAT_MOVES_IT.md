# What moves Coverage@K? (ablations look flat)

**Updated:** 2026-08-15. Evidence from `coverage_curves.py`, Table 8 / D.1 / D.2 COVID CSVs, and C.3 horizon jobs **7627191–95**.

## Definition (code)

Paper **Cov@100** = `coverage_obs_e2` at **K=100** (not the legacy fractional `coverage` / `eps_frac=0.02`).

From `benchmarks/metrics/sequences.py` → `coverage_at_e`:

\[
\mathrm{Cov}@K,\varepsilon
= \frac{|\{t \in \mathrm{obs\ leaves}:\ \exists g\in\mathrm{gen\ pool},\ d_H(t,g)\le \varepsilon\}|}{|\mathrm{obs\ leaves}|}
\]

- **ε** = absolute Hamming (AA edits). Table 5 / C.2 / C.3 paste use **ε=2**.
- **K** = number of generated trees per root; gen leaves are pooled over the first K trees (`coverage_curves.score_pool`).
- **n roots** = typically **5** locked Table-5 roots (COVID Brazil groups 4/5/7/9/10; same for T7/T8/D.*).
- N=16 leaves/tree. So one Cov@100 number averages ≤80 observed leaves × 5 roots.

Legacy `coverage` (fractional ε·L) saturates at **1.0** for Spike L=1280 at ε=0.02 and is **not** the paper column.

## Why T8 / D.1 / D.2 all hit ~0.775

On the **same 5 COVID roots**, NeutralBD, pLM±fit, AR, TreeSBM full, and every Table-8 ablation (`no_bridge`, `no_entropy`, `no_lit_mask`, …) report **`coverage_obs_e2 ≈ 0.775`** at K=100.

| Evidence | Value |
|---|---|
| T8 wide @ K=100 e2 | full / no_bridge / no_entropy / no_bl / no_tree_ctx / no_seq_branch / no_internal all **0.775**; no_lit_mask **0.788** (noise) |
| D.1 / D.2 NeutralBD + ESM2±fit + T8 full/no_bridge | **0.775** |
| COVID Table5 eabs K=10→100 | e1/e2 **flat at 0.775** already by K=10 |

**Why flat:**

1. **Easy leaves dominate.** Most observed leaves on these roots are within 0–2 AA of *some* gen leaf even for NeutralBD (mean_min_edit ≈ 1.2–1.4). Methods that mutate differently still cover the same easy set.
2. **Hard root (group 9 / NODE_0000000)** sits at mean edit ≈5–6 → uncovered at ε=2 for everyone → pulls the mean to the same plateau (~0.775 = majority covered, one hard root uncovered).
3. **Ablations change mutation *identity*, not neighborhood coverage.** Mut recall / aa|hit / Cons / min_edit move; ε=2 set-cover does not.
4. **K=10→100 barely helps at ε=1–2** (already saturated); small gains appear only at ε=0 / ε=3–5 for some ablations (`table8_ablations.md`).

So Cov@100 on COVID Brazil is a **weak discriminator** for reference-process / bridge / architecture ablations. Prefer Mut / aa|hit / Cons / antigenic / min_edit.

## What actually moves coverage

| Lever | Direction | Evidence |
|---|---|---|
| **ε (absolute)** | Strong | COVID T8: e0≈0.60 → e1/e2≈0.775 → e3≈0.79–0.81 → e5≈0.96 → e8≈0.99 |
| **Virus / diversity** | Strong | H3N2 e2≈0.66–0.73; HIV geo e2≈0.19; HIV temporal e2=**0** (mean edit ≈180) |
| **Horizon (genetic H)** | Strong | C.3: COVID med Cov@100=0.938 vs long=0.667; H3 long=**0**; H1 short=1.0 vs long≈0.4–0.5 |
| **K** | Weak at mid-ε on COVID | K=10 vs 100: e2 unchanged; tiny e0/e3/e5 bumps for some ablations. **Cov@500** queued to test high-K (jobs below) |
| **Method / bridge / entropy / lit mask / pLM±fit** | Negligible on COVID e2 | All ≈0.775; mut/aa|hit **do** move (no_bridge mut 0.008 vs full 0.080) |
| **n roots / root lock** | Structural | Same 5 roots → same hard leaf set → same plateau across ablations |

## Practical guidance for the paper

- Report Cov@ε=2 for protocol consistency, but **do not** claim ablation wins from Cov@100 on COVID.
- Show **ε-sweeps** (or e0/e5) and **sequence metrics** when comparing bridge / pLM / NeutralBD.
- Use **C.3 horizon** and **cross-virus** tables when arguing coverage sensitivity.
- Empty C.3 buckets stay empty under locked roots (not a coverage bug).

## Related artifacts

- Def: `benchmarks/coverage_curves.py`, `benchmarks/metrics/sequences.py`
- Flat ablations: `table8_ablations.md`, `coverage_curves_covid_N16_table8_*_K100.csv`
- Horizon: `table_c3_horizon.md` (jobs 7627191–95)
- Cov@500 (C.2): pasted in `table_c2_future_lineage.tex` (COVID TreeSBM still —; resume **7628571**)
