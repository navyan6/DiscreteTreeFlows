# Appendix E.1 / E.2 — Ablations mapped from repo Table 8 (Day 2; E.1 pull 2026-08-17; E.2 queued)

**Updated:** 2026-08-17 (E.1 terminal-only / w/o terminal filled; no_doob still running; E.2 mut/stop-head trains submitted).  
**Source:** [`table8_ablations.md`](table8_ablations.md) + [`PASTE_TABLES_2026-08-14.md`](PASTE_TABLES_2026-08-14.md) §1 (COVID Brazil geo, ckpt `covid_v5_mutrec`, mrs=0.5, n_trees=20).  
**Rule:** map existing T8 cells; do **not** invent factorial rows. Unmapped / unfinished variants stay **—** with blocker.  
**Job IDs:** [`table_e1_job_ids.json`](table_e1_job_ids.json) · [`table_e2_job_ids.json`](table_e2_job_ids.json)

**Additive (on top of paper cols):** Cons↑, Antigenic↑ (PMC lit), AA|hit↑. Also keep Mut / pLM NLL where available.

Cov@100 = `coverage_obs_e2` @ K=100 (**0.775** for almost every ablation on the same 5 roots — weak discriminator). Prefer Mut / aa|hit / min_edit.

Tree-KL† = sim_neutral ≈ ln2 (saturated). Prefer Split-KL from T8b when present.

---

## E.1 Bridge matching ablations

**Virus / protocol:** COVID Brazil geo-test · `covid_v5_mutrec` recipe · mrs=0.5 · n_trees=20 · max_roots=5 (same as T8 / other appendix E/D).

**Paper columns:** Bridge / Doob / Terminal flags + Tree-KL↓, Cov@100↑, Mut. recall↑ (+ Cons, Antigenic, AA|hit)

| Variant | Map / flag | Tree-KL | Cov@100 | Mut | Cons | Antigenic | AA\|hit | pLM NLL |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Reference process only | `no_bridge` (gen `--ablate-bridge`) | 0.693† | 0.775 | 0.008 | 0.876 | 0.004 | 0.086 | 0.4465 |
| Terminal-only generator | train `--ablate-terminal-only` **7628589** | — | **0.775** | **0.089** | **0.797** | **0.009** | **0.231** | **0.4609** |
| Bridge w/o Doob form | train `--ablate-doob` **7628590** RUNNING | — | — | — | — | — | — | — |
| TreeSBM w/o terminal consistency | train `--ablate-terminal-consistency` **7628591** | — | **0.775** | **0.062** | **0.794** | **0.008** | **0.201** | **0.4599** |
| TreeSBM full | `full` | 0.693† | 0.775 | 0.080 | 0.839 | 0.005 | 0.375 | 0.4572 |

### What each ablation means

| Variant | Code effect | Why retrain (not gen-time) |
|---|---|---|
| Reference only | Gen: force `log R_θ = log R0` | Already T8 `--ablate-bridge` on full ckpt |
| Terminal-only | Train KL target forced to t→1 pure CE on x1 (no bridge mixture at t\<1) | Changes loss target → new `c_θ` |
| Bridge w/o Doob | Train KL(softmax(R0) \|\| R_θ); skip Doob h-transform / ignore x1 | Changes loss target → new `c_θ` |
| w/o terminal consistency | Train `lambda_cons=0` (drop L_cons anchoring) | Changes loss weights → new `c_θ` |
| Full | Standard Doob bridge + L_mut + L_cons | Existing `covid_v5_mutrec` |

### Betty jobs (submitted 2026-08-15)

| Variant | Train | Enrich (afterok) | Cov@100 (afterok) | Ckpt dir |
|---|---:|---:|---:|---|
| Terminal-only | **7628589** | **7628592** | **7628593** | `checkpoints/covid_e1_terminal_only` |
| w/o Doob | **7628590** | **7628594** | **7628595** | `checkpoints/covid_e1_no_doob` |
| w/o terminal | **7628591** | **7628596** | **7628597** | `checkpoints/covid_e1_no_terminal` |

Script: `scripts/slurm_e1_bridge_ablations.sh` (`ABLATE=terminal_only|no_doob|no_terminal`, `MODE=train|enrich|coverage|baselines`).

**Timelimit:** train jobs **7628589–91** bumped via `scontrol` to **72h** (`3-00:00:00`); IDs unchanged; afterok enrich/cov **7628592–97** intact.  
**Pull 2026-08-17:** terminal-only train/en/cov **COMPLETED**; w/o terminal train/en/cov **COMPLETED**; w/o Doob train **7628590 RUNNING** (~1-21h / 72h) so enrich/cov **7628594/95** still Dependency. Tree-KL baselines not queued (saturates ln2). Do **not** invent no_doob numbers.

**Related T8 cells (not paper E.1 rows; optional caption):**

| T8 ablation | Cov@100 | Mut | Cons | Ant | AA\|hit | Note |
|---|---:|---:|---:|---:|---:|---|
| w/o seq-dependent branching | 0.775 | 0.038 | 0.818 | 0.011 | 0.242 | closer to “pLM + branching off”; see Day 4 [`table_d2_ref_components.md`](table_d2_ref_components.md) |
| w/o lit/PMC mask | 0.788 | 0.069 | 0.839 | —¶ | 0.341 | gen-time mask drop; enrich ant null |

¶ Antigenic — in no_lit_mask enrich (mask dropped for scoring).

Jobs: no_bridge enrich **7597319** / cov **7585972** / KL **7585977**; full enrich **7585970** / cov **7585971** / KL **7585976**.

---

## E.2 Architecture ablations

**Paper columns:** Tree encoder / Branch / Mutation / Stop heads + Tree-KL↓, Cov@100↑ (+ Cons, Antigenic, AA|hit, Mut, pLM NLL)

| Variant | Map from T8 | Tree-KL | Cov@100 | Mut | Cons | Antigenic | AA\|hit | pLM NLL |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Node-independent rates | `no_tree_ctx` | —‡ | 0.775 | 0.084 | 0.845 | 0.012 | 0.264 | 0.4562 |
| No branch-length head | `no_bl_eval` (gen-time) | 0.693† | 0.775 | 0.080 | 0.839 | 0.005 | 0.375 | 0.4572 |
| No mutation head | train `--ablate-mut-head` **7642530** (+en **7642532** / cov **7642533**) | — | — | — | — | — | — | — |
| No stop/termination head | train `--ablate-stop-head` **7642531** (+en **7642534** / cov **7642535**) | — | — | — | — | — | — | — |
| Full TreeSBM | `full` | 0.693† | 0.775 | 0.080 | 0.839 | 0.005 | 0.375 | 0.4572 |

† Tree-KL saturated. ‡ **Blocker:** KL job not run for `no_tree_ctx` (would be ln2 anyway). Split-KL — same.

**E.2 train flags (new 2026-08-17):** `--ablate-mut-head` forces `log R_θ = log R0` and `λ_mut=0` (branching/BL/stop still trained); `--ablate-stop-head` uses constant `p_stop=0.5` and `λ_stop=0`. Ckpt dirs `checkpoints/covid_e2_{no_mut_head,no_stop_head}/` with `--resume`. 72h train wall. Coverage + pLM NLL are additive evals (Tree-KL not queued). Alg. 4 does not sample the stop head; stop ablation is train-loss / encoder only. IDs: [`table_e2_job_ids.json`](table_e2_job_ids.json). Leave metric cells **—** until afterok enrich/cov land.

**Additional T8 architecture-adjacent (caption only, not paper E.2 rows):**

| T8 ablation | Cov@100 | Mut | Cons | Ant | AA\|hit | pLM NLL | Tree-KL / Split-KL |
|---|---:|---:|---:|---:|---:|---:|---|
| w/o per-site entropy | 0.775 | 0.085 | 0.835 | 0.005 | 0.312 | 0.4573 | 0.693 / 57.61 |
| w/o internal-node seqs | 0.775 | 0.090 | 0.834 | 0.006 | 0.270 | 0.4565 | — / — |
| λ_br=0 retrain **7539981** | — | — | — | — | — | — | **TIMEOUT** — do not paste |

Jobs: no_tree_ctx enrich **7597320** / cov **7597326**; no_bl_eval enrich **7597322** / cov **7585974** / KL **7585979**.

---

## Caption guidance

- E.1/E.2 are **mapped from COVID Table 8**, not a fresh factorial.  
- E.1 terminal-only / w/o Doob / w/o terminal are **train-time** loss ablations (not gen-time `--ablate-*` on `covid_v5_mutrec`).  
- E.2 no-mut / no-stop heads are **train-time** (`--ablate-mut-head` / `--ablate-stop-head`); ckpt config is applied at load so enrich/coverage need no extra gen flags.  
- Cov@100 flatness: same 5 Brazil roots; group 9 hard (mean edit≈5.6).  
- Fitness ablation is **not** a T8 row — mapped from D.1 / paper Table 6 on Day 4 (D.2 → [`table_d2_ref_components.md`](table_d2_ref_components.md)).
