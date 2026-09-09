# Paste-ready tables (2026-08-14)

Last refresh: **2026-08-25** (COVID TreeSBM Cov@500 = **0.775** from K500 CSV; COVID pLM EVEscape Δ = **−0.424** from `eval_table5_covid_plm_evescape.json`). Do not invent numbers. Dashes = unexplained below or jobs pending.

---

## 0. OAS / Ab Track C status (user-facing)

**Prior break.** FastTree bootstrap support glued onto Newick node names → `*_bl.json` key mismatch → `Some nodes are not reachable by root`. Fixed by `scripts/ab_postprocess_clones_for_treesbm.py` (611/611 OK).

**Gap-fill break (7618018).** `anc_aa` had root+tips only; missing internals filled with `"-"*L` → bridge → `val=nan` → no `best.pt`.

**Betty now (2026-08-15 Day 1):**

| Job | Role | State |
|---:|---|---|
| **7627289** | OAS 1M TreeSBM train (`ab_oas_1m_v1`) | **FAILED** exit 1 — but **`best.pt` written** (best val=**64.4053**, early-stop @ep122). Crash *after* ckpt: post-train eval missing PLM cache |
| **7627915** | OAS PLM precompute | **COMPLETED** — train 397 / val 44 / test 170 complete groups cached |
| **7627919** | Rod.82 eval as `treesbm_ab_oas` | **COMPLETED** 01:03:40 → metrics in `table6_ab_track_c.md` |
| 7627290 / 7627916 | prior Rod.82 attempts | CANCELLED / FAILED (`aa_indices`) |

Do **not** mix Track A OAS with Rodriguez-82 Track B/C metrics. Paper Table 5 / Track C paste includes **treesbm_ab_oas** (below §3b).

---

## 1. Paper “Table 7” = Core TreeSBM ablations (repo Table 8) — COVID Brazil geo

ckpt=`covid_v5_mutrec` · mrs=0.5 · n_trees=20 · Cov@100 = `coverage_obs_e2` · Antigenic = PMC lit hotspot frac.

### Paste-ready

| Variant | Tree-KL | Cons | AA\|hit | Antigenic (PMC) | Cov@100 | Mut recall | pLM NLL |
|---|---:|---:|---:|---:|---:|---:|---:|
| w/o bridge | 0.6931† | 0.8758 | 0.0858 | 0.0043 | 0.775 | 0.0082 | 0.4465 |
| w/o fitness weighting in reference | —‡ | —‡ | —‡ | —‡ | —‡ | —‡ | —‡ |
| w/o seq-dependent branching | —† | 0.8181 | 0.2418 | 0.0114 | 0.775 | 0.0375 | 0.4540 |
| w/o BL head | 0.6931† | 0.8390 | 0.3749 | 0.0051 | 0.775 | 0.0795 | 0.4572 |
| w/o tree-context | —† | 0.8448 | 0.2635 | 0.0120 | 0.775 | 0.0841 | 0.4562 |
| w/o per-site entropy | 0.6931† | 0.8346 | 0.3121 | 0.0053 | 0.775 | 0.0845 | 0.4573 |
| full | 0.6931† | 0.8389 | 0.3753 | 0.0051 | 0.775 | 0.0796 | 0.4572 |

† **Tree-KL:** empirical track is NaN (Dirac JS vacuous without matched clade sets in older protocol; current empirical KL still NaN in CSV). Values shown are **sim_neutral Tree-JS ≈ ln2 = 0.693147** from T8b jobs **7585976–80** — saturated for every ablation, so **not discriminative**. Prefer **Split-KL** if you need a topology distribution number:

| Ablation | Split-KL (sim_neutral N=16) | job |
|---|---:|---:|
| no_bridge | 57.44 | 7585977 |
| no_entropy | 57.61 | 7585978 |
| no_bl_eval | 57.54 | 7585979 |
| no_lit_mask | 57.54 | 7585980 |
| full | 57.47 | 7585976 |
| no_tree_ctx / no_seq_branch | — | Cov/KL jobs not run (only enrich + late Cov **7597326/27**) |

‡ **Fitness row:** not a T8 gen-time ablation. T8 “full” already has `fitness_beta` unset (β=0). Closest mapped numbers live in **Table 6 / repo Table 7** (ESM-2-650M ±fit, or TreeSBM β=1). Do not paste a fake T8 row.

**Cov@100 flat 0.775:** same 5 Brazil roots (groups 4,5,7,9,10). Group 9 is hard (mean edit≈5.6 → Cov@e2=0 on that root); groups 5/7/10 are near-identical leaves → overall mean ≈ 0.775 for almost every method. Ablations differ in **mut / aa|hit / min_edit**, not mid-ε coverage. New Cov@100 for tree-ctx / seq-branch / internal: all e2=**0.775** (jobs **7597326–28** COMPLETED).

**pLM NLL:** filled from re-enrich **7597319–25** (+ prior 7585969/70). Values ≈0.45–0.46 except no_bridge **0.4465**.

Sources: `checkpoints/table8_*_mrs0.5.json` (nested `summary`), `coverage_curves_covid_N16_table8_*_K100.csv`.

---

## 2. Paper “Table 6” = Reference process (repo Table 7) / appendix D.1

### Paste-ready

| R0 / method | Tree-KL | Split-KL | Cov@100 | Mut recall | AA\|hit | Antigenic (PMC) | Cons | pLM NLL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| JTT/WAG/LG (no fit) | — | — | 0.775 | 0.0223 | 0.0825 | 0.0066 | 0.8213 | **0.4705** |
| ESM-2-650M (no fit) | — | — | 0.775 | 0.0180 | 0.1083 | 0.0059 | 0.8931 | **0.4415** |
| ESM-2-650M (fit) | — | — | 0.775 | 0.0011 | 0.0983 | 0.0018 | 0.9884 | **0.4405** |
| ESM-C (no fit) | — | — | —¶ | 0.0347 | 0.0796 | 0.0055 | 0.6923 | **0.4605** |
| ESM-C (fit) | — | — | —¶ | 0.0356 | 0.0853 | 0.0053 | 0.7134 | **0.4267** |
| ProGen2 (fit) | — | — | —¶ | 0.1126 | 0.3978 | 0.0048 | 0.9855 | **0.4546** |
| TreeSBM default (β=1) | 0.6931† | 57.52† | 0.775 | 0.0378 | 0.2896 | 0.0061 | 0.9412 | **0.4464** |
| Evo2 | — | — | — | — | — | — | — | — | **Blocker:** Evo2 backend not wired |

† Tree-KL/Split-KL only for TreeSBM from all-method COVID baselines **7575017** sim_neutral N=16 — **not** per-R0. Tree-KL again saturates at ln2. Per-backend KL would need `MODE=baselines` per ROW (not submitted; expensive and Tree-KL still saturates).

¶ **Cov@100 Day 5:** jobs **7628260** / **7628261** / **7628262** — Betty SSH **port 22 timeout** (~13:31 ET Day 5 retry). Leave Cov cells **—**; do not invent; do not resubmit until sacct confirms FAILED. Enrich already filled. Pull when Betty reachable.

**ESM-C / ProGen2 (2026-08-15):**
- **ESM-C:** HF `Synthyra/ESMplusplus_small`. Jobs **7627569** (nofit) / **7627570** (fit) — full n=20.
- **ProGen2:** enijkamp `progen2-small`. Job **7627918** — **COMPLETED n=20** (supersedes partial **7627571** n=3). Artifact: `checkpoints/table7_progen2_fit_mrs0.5.json`.
- Full detail: `table7_refproc.md`.

**pLM NLL:** all T7 rows filled from enrich JSONs (`checkpoints/table7_*_mrs0.5.json`). Jobs **7626622–25**, **7627569–70**, **7627918**.

---

## 3. Table C.1 — Full tree generation across datasets

**Full Day-2 paste:** [`table_c1_full_tree_gen.md`](table_c1_full_tree_gen.md).

Columns: Tree-KL, Split-KL, RF, Quartet, Branch W1, pLM NLL, aa|hit, antigenic, cons.

**Tree-KL honesty:** empirical Tree-KL/Split-KL = **NaN** in CSV. sim_neutral Tree-KL ≈ **ln2** for all methods (vacuous saturation). Prefer RF / Quartet / W1 on empirical; Split-KL on sim if needed.

### Paste-ready (filled cells only)

#### Influenza HA (H3N2) — empirical N=16 (`results_baselines.csv` mean)

| Method | Tree-KL | Split-KL | RF | Quartet | Branch W1 | pLM NLL | aa\|hit | Antigenic (flu lit) | Cons |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Neutral CTMC | 0.693† | 43.0‡ | 0.980 | 0.714 | 7.77e-4 | — | — | — | — |
| Autoregressive (adapted) | 0.693† | 62.9‡ | 0.986 | 0.717 | 9.21e-4 | — | — | — | — |
| TreeSBM | 0.693† | 57.4‡ | 0.983 | 0.721 | 9.25e-4 | — | **0.439**¶ | **0.006**¶ | **0.918**¶ |

†‡ sim_neutral N=16 means. ¶ `eval_enrichment_h3n2_v3_lit_hotspot_mrs0.5_fullmetrics.json`. Neutral/AR seq — (no enrich).

#### SARS-CoV-2 Spike — empirical N=16 (`results_baselines_covid.csv` job **7575017**)

| Method | Tree-KL | Split-KL | RF | Quartet | Branch W1 | pLM NLL | aa\|hit | Antigenic | Cons |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Neutral CTMC | 0.693† | 43.0‡ | 0.980 | 0.672 | 2.29e-5 | — | — | — | — |
| Autoregressive (adapted) | 0.693† | 62.9‡ | 0.990 | 0.668 | 5.62e-5 | — | 0.094§ | 0.017§ | 0.993§ |
| TreeSBM | 0.693† | 57.5‡ | 0.982 | 0.669 | 5.45e-5 | 0.457¶ | 0.375¶ | 0.005¶ | 0.839¶ |

§ `eval_table5_covid_baselines.json`. ¶ T8 full enrich. Neutral seq — (blocker).

#### Influenza HA H1N1 — geo N=16 (replaces HIV Env in C.1 `.tex`)

**Best TreeSBM:** `h1n1_v2_lit_hotspot` · geo `data/h1n1` (Table 5 / C.3). **Not** leafholdout / temporal.

| Method | Tree-KL | Split-KL | RF | Quartet | Branch W1 | pLM NLL |
|---|---|---|---|---|---|---|
| Neutral CTMC | -- | -- | -- | -- | -- | -- |
| Autoregressive tree-edit | -- | -- | -- | -- | -- | -- |
| TreeSBM | -- | -- | -- | -- | -- | -- |

**Blocker / jobs (do not invent):** no prior `results_baselines_h1n1.csv`; T5 JSON **7539950** has no `plm_nll`. Queued **7646375** topo · **7646376** Neutral/AR NLL · **7646377** TreeSBM NLL. IDs: [`table_c1_h1n1_job_ids.json`](table_c1_h1n1_job_ids.json). HIV Env numbers remain in `hiv_{geo,temporal}/` but are **out of C.1 LaTeX**.

#### Antibody / Protein families

| Row | Status |
|---|---|
| Ab Track C distributional | **DONE** — §3b |
| Ab C.1 RF/Quartet/W1 | **FILLED** free-topo **7635753** |
| Protein families | **—** no PFAM jobs |

---

## 3b. Paper “Table 5” Ab / Track C — Rodriguez 82

Writeup: [`table5_ab_affinity_maturation.md`](table5_ab_affinity_maturation.md). Full grid: [`table6_ab_track_c.md`](table6_ab_track_c.md).

### Paste-ready (matches paper screenshot columns; verified numbers)

| Method | CDR mut. recall ↑ | SHM load error ↓ | Terminal diversity error ↓ | Coverage@100 ↑ |
|---|---:|---:|---:|---:|
| Neutral SHM model | 0.688 | 0.054 | 17.435 | 0.059 |
| Autoregressive tree-edit model | 0.000 | 0.117 | 19.142 | 0.081 |
| CTMC model (CoSiNE) | 0.658 | 0.037 | 3.019 | 0.076 |
| TreeSBM (pathogen `best.pt`) | 0.278 | 0.055 | 5.104 | 0.084 |

- **Coverage@100** = absolute Hamming **ε=5**, K=100 (not viral ε=2).
- **Neutral ≠ Thrifty.** Data = Rodriguez 82 frozen families; TreeSBM row = OOD pathogen ckpt.
- Screenshot typos: CoSiNE CDR 0.678→**0.658**; AR Cov 0.0417→**0.081**; TreeSBM term-div 9.088→**5.104**.

---

## 3c. Appendix C.3 — Forecast horizon

Full paste: [`table_c3_horizon.md`](table_c3_horizon.md). Jobs **7627191–95** COMPLETED. Tree-KL column **—** (deferred). Empty buckets locked empty. Additive Cons / AA|hit in CSV; **Antigenic filled** from cache score **7628578** (PMC / flu lit / HIV V1–V5). Locked-root explanation in protocol. Cov flatness: [`COVERAGE_WHAT_MOVES_IT.md`](COVERAGE_WHAT_MOVES_IT.md).

| Dataset | Horizon | Method | Cov@100 | Mut | Clade | Min dist | Tree-KL | Cons | AA\|hit |
|---|---|---|---:|---:|---:|---:|---|---:|---:|
| SARS-CoV-2 | medium | Neutral/AR/TreeSBM | 0.938 | ~0 | 0 | 0.69 | — | 1.000 | —/—/— |
| SARS-CoV-2 | long | Neutral/AR/TreeSBM | 0.667 | ~0 | 0 | ~1.9 | — | 1.000 | sparse |
| H1N1 | short | Neutral/AR/TreeSBM | 1.000 | 0.045/0.067/0.212 | 0.75/0.25/0.50 | 0.44/0.46/0.31 | — | 1.000 | 1.000 |
| H1N1 | long | Neutral/AR/TreeSBM | 0.406/0.406/0.500 | 0.013/0.005/0.141 | 0/0/0.10 | 2.91/2.91/2.41 | — | 1.000 | 1.000 |
| H3N2 | medium | Neutral/AR/TreeSBM | 0.828/0.859/0.859 | 0/0.014/0.099 | 0/0.25/0 | 1.27/1.22/1.14 | — | 1.000 | —/1/1 |
| H3N2 | long | Neutral/AR/TreeSBM | 0.000 | 0/0/0.131 | 0 | 8.25/8.25/7.19 | — | 1.000 | —/—/1 |
| HIV geo | short | Neutral/TreeSBM | 0.313 | 0.001 | 0.067/0 | 88.6 | — | 1.000 | 0.523/0.740 |
| HIV geo | long | Neutral/TreeSBM | 0.000 | 0.023/0.050 | 0.097 | 26.8/24.7 | — | 0.982/0.994 | 0.576/0.695 |
| HIV temporal | medium | Neutral/TreeSBM | 0.000 | 0.016/0.020 | 0/0.113 | 177.6/177.1 | — | 1.000 | 0.694/0.701 |
| HIV temporal | long | Neutral/TreeSBM | 0.000 | 0.023/0.037 | 0.067/0.098 | 186.5/184.3 | — | 0.973/0.994 | 0.212/0.400 |

---

## 3d. Appendix J.1 — Linear-probe (FINAL)

Confirmed in [`PAPER_TABLES.md`](../../PAPER_TABLES.md) (H3N2 job **7434252**, COVID job **7431998**). No new Day-1 compute. Primary paste already FINAL — paper table is H3N2-shaped; COVID companion optional.

---

## 3e. Appendix C.2 — Future-lineage (Day 2 + Cov@500 pull 2026-08-16)

LaTeX: [`table_c2_future_lineage.tex`](table_c2_future_lineage.tex). Full notes: [`table_c2_future_lineage.md`](table_c2_future_lineage.md). **Naming:** user “c.3” paste = **C.2** (`tab:app_future_lineage_full`); horizon = C.3.

**Paper rows (pLM + TreeSBM only):**

| Dataset | Method | Cov@100 | Cov@500 | Mut | Clade | Min dist |
|---|---|---:|---:|---:|---:|---:|
| SARS-CoV-2 Spike | pLM prior | 0.775 | **0.775** | 0.048 | 0.000 | 1.425 |
| SARS-CoV-2 Spike | TreeSBM (`covid_v5_mutrec`) | 0.788 | **0.775** | 0.080 | 0.000 | 1.388 |
| Influenza HA H3N2 | pLM prior | 0.675 | **0.750** | 0.004 | 0.000 | 2.650 |
| Influenza HA H3N2 | TreeSBM (`h3n2_v3_lit_hotspot`) | 0.725 | **0.788** | 0.100 | 0.100 | 2.238 |
| HIV Env | pLM prior | **0.188** | **0.188** | **0.015** | **0.064** | **63.9** |
| HIV Env | TreeSBM (**geo**) | 0.188 | **0.200** | 0.026 | 0.061 | 63.1 |
| Antibody lineages | pLM prior | 0.018 | **0.018** | CDR **0.297** | **0.039** | **15.9** |
| Antibody lineages | TreeSBM (`treesbm_ab_oas`) | **0.027** | **0.027** | CDR **0.150** | **0.018** | **14.3** |

Cov@500: all paper rows filled (COVID TreeSBM **0.775** = plateau at mid-ε). Coverage flatness: [`COVERAGE_WHAT_MOVES_IT.md`](COVERAGE_WHAT_MOVES_IT.md).

**Additive context (AR / Neutral / H1 / temporal — not in bare tex):**

| Dataset | Method | Cov@100 (e2) | Mut | Clade | Min dist | Cons | Antigenic | AA\|hit |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| COVID | pLM / AR / TreeSBM | 0.775 / 0.775 / 0.788 | 0.048 / 0.006 / 0.080 | 0 / 0 / 0 | 1.425 / 1.425 / 1.388 | 0.985 / 0.993 / 0.839 | 0 / 0.017 / 0.005 | 0.017 / 0.094 / 0.375 |
| H1N1 | pLM / AR / TreeSBM | 0.762 / 0.762 / 0.762 | 0.071 / 0.006 / 0.241 | 0 / 0 / 0.175 | 1.462 / 1.462 / 1.388 | 0.945 / 0.998 / 0.952 | 0.115 / 0 / 0.163 | 0.136 / 0.008 / 0.618 |
| H3N2 | Neut / pLM / AR / TreeSBM | 0.662 / 0.675 / 0.675 / 0.725 | 0.020† / 0.004† / 0.079† / **0.100** | 0 / 0 / 0.30 / 0.10 | 2.65 / 2.65 / 2.40 / 2.24 | 1† / 1† / 1† / **0.918** | — / — / — / **0.006** | — / — / — / **0.439** |
| HIV geo | Neut / AR / TreeSBM | 0.188 | 0.002 / 0.025 / 0.026† | 0.079 / 0.087 / 0.061 | ~63 | 0.786 / 0.984 / 0.893 | 0.121 / 0.160 / 0.134 | 0.146 / 0.125 / 0.337 |
| HIV temporal | Neut / AR / TreeSBM | **0.000** | 0.012 / 0.024 / 0.034† | 0.060 / 0.104 / 0.040 | ~182 | 0.728 / 0.974 / 0.870 | 0.126 / 0.148 / 0.099 | 0.071 / 0.184 / 0.359 |
| Ab Rod.82 | pLM / **treesbm_ab_oas** | 0.018 / **0.024** | CDR 0.297 / **0.150** | — | — | — | — | — |

† H3 Neut/pLM/AR mut+cons from coverage CSV (no enrich); TreeSBM additives from enrich. HIV mut† from coverage; Cons/Ant/aa|hit from enrich.

---

## 3f. Appendix E.1 / E.2 — mapped from Table 8 (Day 2)

Full paste: [`table_e1_e2_ablations.md`](table_e1_e2_ablations.md).

| E.1 / E.2 row | Map | Cov@100 | Mut | Cons | Ant | AA\|hit | pLM NLL |
|---|---|---:|---:|---:|---:|---:|---:|
| E.1 Reference only | T8 no_bridge | 0.775 | 0.008 | 0.876 | 0.004 | 0.086 | 0.4465 |
| E.1 Terminal-only | `--ablate-terminal-only` **7628589/92/93** | **0.775** | **0.089** | **0.797** | **0.009** | **0.231** | **0.4609** |
| E.1 w/o Doob | `--ablate-doob` **7628590** RUNNING | — | — | — | — | — | — |
| E.1 w/o terminal | `--ablate-terminal-consistency` **7628591/96/97** | **0.775** | **0.062** | **0.794** | **0.008** | **0.201** | **0.4599** |
| E.1 TreeSBM full | T8 full | 0.775 | 0.080 | 0.839 | 0.005 | 0.375 | 0.4572 |
| E.2 Node-independent | T8 no_tree_ctx | 0.775 | 0.084 | 0.845 | 0.012 | 0.264 | 0.4562 |
| E.2 No BL head | T8 no_bl_eval | 0.775 | 0.080 | 0.839 | 0.005 | 0.375 | 0.4572 |
| E.2 No mut head | `--ablate-mut-head` train **7642530** (+en **7642532** / cov **7642533**) | — | — | — | — | — | — |
| E.2 No stop head | `--ablate-stop-head` train **7642531** (+en **7642534** / cov **7642535**) | — | — | — | — | — | — |
| E.2 Full | T8 full | 0.775 | 0.080 | 0.839 | 0.005 | 0.375 | 0.4572 |

---

## 3g. Appendix D.2 — Reference-process components (Day 4)

Full paste: [`table_d2_ref_components.md`](table_d2_ref_components.md). Mapped from NeutralBD + T7 ESM2±fit + T8. **New jobs:** Branch W1 ESM2 nofit/fit **7628576/77**; Neutral ant/aa|hit via cache score **7628578**. Tree-KL not pasted for ESM2 (saturated). Caption as overlapping map, not dedicated factorial.

| Variant | Cov@100 | Tree-KL | Branch W1 | Cons | Antigenic | AA\|hit |
|---|---:|---:|---:|---:|---:|---:|
| Neutral birth-death | 0.775 | 0.693† | 2.29e-5 | 1.000 | **0.012** | **0.828** |
| pLM mutation only (ESM2 nofit) | 0.775 | — | **5.20e-5** | 0.893 | 0.006 | 0.108 |
| pLM + fitness (ESM2 fit) | 0.775 | — | **5.32e-5** | 0.988 | 0.002 | 0.098 |
| pLM + branching (T8 full β=0) | 0.775 | 0.693† | 5.26e-5 | 0.839 | 0.005 | 0.375 |
| pLM + fit + branching (T8 no_bridge‡) | 0.775 | 0.693† | 5.49e-5 | 0.876 | 0.004 | 0.086 |
| TreeSBM full bridge (T8 full) | 0.775 | 0.693† | 5.26e-5 | 0.839 | 0.005 | 0.375 |

† Saturated ln2. ‡ Careful caption: no_bridge = R0 / ablate-bridge. Neutral ant/aa from cache **7628578**. ESM2 W1 from **7628576/77** COMPLETED. Main Table 8 panels = **—**.

---

## 3h. Appendix E.3 — Sampling sensitivity (2026-08-16 grid)

Full paste: [`table_e3_sampling.md`](table_e3_sampling.md). COVID one-factor grid: K∈{50,100,500}, T∈{0.5,1.0,1.5}, max_leaves∈{100,400,800}. Jobs **7639276–79** (temp 0.5/1.5 + maxl 100/800); K500 via C.2 **7635911**.

| Setting | Cov | Diversity | Tree-KL | Runtime | Cons / Mut |
|---|---:|---|---|---:|---|
| H3 K=50 @e5 | 0.812 | — | — | 15668 s / 5 roots | cons≈1; mut=0.252 |
| H3 K=100 @e5 | 0.812 | — | — | same | mut=0.156 |
| COVID K=50 @e2 | 0.775 | — | — | 21421 s / 5 roots | — |
| COVID K=100 @e2 | 0.788 | — | — | same | — |
| COVID K=500 | 0.775 (e2) | — | — | — | K500 CSV; mid-ε plateau |
| T=0.5 | **0.775** | — | — | 20498 s | cons=1; mut=0.001 · **7639276** |
| T=1.0 | **0.788** | — | — | 21421 s | Table5 treesbm |
| T=1.5 | **0.775** | — | — | 20639 s | cons=1; mut=0.001 · **7639277** |
| max_leaves=100 | **0.775** | — | — | 24655 s | mut=0.039 · **7639278** |
| max_leaves=400 | **0.788** | — | — | 21421 s | production default / Table5 |
| max_leaves=800 | **0.775** | — | — | 44789 s | mut=0.042 · **7639279** |

---

## 3i. Appendix I.1 — Compute cost (Day 5)

Full paste: [`table_i1_compute.md`](table_i1_compute.md).

| Method | Params | Time/100 (5-root pool) | GPU mem |
|---|---|---|---|
| pLM prior | ESM-2 8M | ~15 s COVID | — |
| AR adapted | — | ~153 s COVID | — |
| NeutralBD | — | ~188 s COVID | — |
| TreeSBM | ~3M + 8M PLM | ~5.95 h COVID / ~4.35 h H3 | — |

---

## 3j. Appendix J.2 / J.3 — Hyperparams & sampling (Day 5)

Full paste: [`table_j2_j3_hyperparams.md`](table_j2_j3_hyperparams.md). Highlights: d_model=128, 4 layers, 8 heads, ESM-2 8M; coverage n_steps=50, enrich n_steps=100, N=16, K_max=100, mrs=0.5, site_temperature=1.0.

---

## 3k. Appendix H.1 / H.2 — Qualitative (Day 5)

Full paste: [`table_h1_h2_qualitative.md`](table_h1_h2_qualitative.md). COVID g40 matched 207-tip; H3 screen #1 g55 / #2 g40. HIV/Ab galleries **—**.

---

## 3l. Paper Table 4 / repo Table 5 — Viral evolution forecasting (EVEscape Δ)

Full paste: [`table5_viral_forecasting.md`](table5_viral_forecasting.md). EVEscape Δ = `model_evescape − gt_evescape` (all scored muts) from `eval_table5_*_baselines.json`.

| Virus | Method | Cov@100 | Mut | EVEscape Δ | Cons | Source |
|---|---|---:|---:|---:|---:|---|
| SARS-CoV-2 | pLM prior | 0.775 | 0.048 | **−0.424** | 0.985 | `eval_table5_covid_plm_evescape.json` (**7644883**) |
| SARS-CoV-2 | AR | 0.775 | 0.006 | **−0.371** | 0.993 | 7539949 |
| SARS-CoV-2 | TreeSBM | 0.788 | 0.080 | **−0.048** | 0.839 | T8 full / 7539980 |
| H1N1 | pLM prior | 0.762 | 0.071 | **−0.354** | 0.945 | 7539950 |

COVID pLM EVEscape was NaN under ESM L=566 (**7539949**); filled from L=1280 re-run.

---

## Job ID cheat sheet (Day 1–5)

| ID | What | State (~13:35 ET 2026-08-15) |
|---:|---|---|
| **7627918** | T7 ProGen2 enrich (n=20) | **COMPLETED** |
| **7628260 / 61 / 62** | T7 Cov@100 ESM-C nofit / fit / ProGen2 | **UNKNOWN** (Day 5 Betty SSH timeout; leave Cov —; no resubmit) |
| **7627569 / 70** | T7 ESM-C enrich | **COMPLETED** |
| **7627915 / 7627919** | OAS PLM precompute / Rod.82 | **COMPLETED** |
| **7627191–95** | C.3 horizon | **COMPLETED** |
| **7626622–25** | T7 pLM NLL re-enrich (JTT/ESM2±/TreeSBM) | **COMPLETED** |
| 7597319–25 / 7597326–28 | T8 pLM NLL re-enrich / remaining Cov | **COMPLETED** |
| **7626610 / 7626609** | HIV geo / temporal baselines | **COMPLETED** |
| **7517208** | Ab thrifty (additive; **not** Neutral SHM) | **COMPLETED** |
| **nohup** | Ab neutral_shm JC69 + ar_tree_edit Rod.82 | **DONE** (login01; sbatch DNS broken) |
| **nohup** | Ab plm_prior Rod.82 | **RUNNING** (~533/1640) |
| **7626612 / 7626613** | HIV Neutral/AR enrich | **COMPLETED** |
| **7626614 / 7626615** | HIV temporal/geo cov | **COMPLETED** |
| 7617937–40 | HIV TreeSBM enrich+cov | **COMPLETED** |
| **7575017** | COVID NeutralBD Tree-KL/W1 (D.2) | **COMPLETED** |
| **7459412** | H3N2 coverage K-sweep (E.3) | **COMPLETED** |
