# Appendix H.1 / H.2 — Qualitative galleries

**Status 2026-08-15 Day 5:** figure **selection** from existing artifacts (no new jobs).

Main-text Figures 2/4/5 already carry qualitative story; appendix galleries optional.

---

## H.1 — Observed vs generated trees (topology)

| Panel | Dataset | Recommended pair | Artifact dir | Rank / note |
|---|---|---|---|---|
| H.1a | SARS-CoV-2 Spike | `group_040_observed.nwk` ↔ `group_040_generated_matched.nwk` (both **207** tips) | `results/covid_tree_viz/` | leaf-count matched; ckpt `covid_v7_pmc_hotspot`, mrs=0.5, job **7437255** |
| H.1b | SARS-CoV-2 Spike (fast) | `group_040_observed_subsample.nwk` ↔ `group_040_generated.nwk` (both **48** tips) | same | jobs **7437188** |
| H.1c | Influenza H3N2 | `group_055_observed_matched.nwk` ↔ `group_055_generated_matched.nwk` | `results/h3n2_tree_viz/` | screen rank **#1** (`SCREEN_RANKING.md`) |
| H.1d | Influenza H3N2 alt | group **40** matched pair | same | screen rank **#2** |

**Do not** compare full observed tip count to unmatched fast gens (README warning).

Screen rankings: `results/covid_tree_viz/SCREEN_RANKING.md`, `results/h3n2_tree_viz/SCREEN_RANKING.md`.

---

## H.2 — Sequence / mutation trajectory callouts

| Panel | Content | Artifact |
|---|---|---|
| H.2a | COVID group **3** Gamma-tip substitution recovery | `results/covid_tree_viz/group_003_gamma_figure_stats.md` (+ `.json` / CSV) |
| H.2b | COVID closest-group screen (group **9** top combo) | `SCREEN_RANKING.md` rank 1; matched gens under `group_009_*` |
| H.2c | Coverage curves (supp figure, not tree gallery) | `benchmarks/results/plots/coverage_curves_h3n2_N16.{png,pdf}` |

---

## Blockers / honesty

| Item | Status |
|---|---|
| Rendered PDF appendix panels | **Author paste** from Newick/stats above into paper figures — repo has trees + rankings, not final illustrator PDFs |
| HIV / Ab qualitative gallery | **—** · **Blocker:** no `*_tree_viz` pack locked for HIV or Rod.82 |
| New generations | **Not started** (Day 5 forbids new jobs) |
