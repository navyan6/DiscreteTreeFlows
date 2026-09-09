# NeurIPS TreeSBM paper review

**PDF:** `TreeSBM__NeurIPS___1_.pdf` (14 pp; pdfTeX 2026-08-18 03:40 UTC). Extracted with pypdf (text quality adequate for all tables/captions).  
**Rule:** no invented numbers. Paper values vs repo artifacts only.  
**Authoritative sources:** `APPENDIX_COMPLETION_PLAN.md`, `PASTE_TABLES_2026-08-14.md`, `TABLE_METADATA.md`, `PAPER_TABLES.md`, `COVERAGE_WHAT_MOVES_IT.md`, `antibody_tracks.md`, and the `table_*.md`/`.tex` companions listed below.  
**Numbering map:** this NeurIPS PDF does **not** use the ICLR appendix IDs. Paper E.1 ≈ repo D.1 / Table 7; paper E.2 ≈ repo D.2; paper F.1 ≈ repo Table 8 / E.2 architecture. Cited below as **Paper T# / Paper App X**.

---

## 1. Executive summary

The PDF is a 6-page main + 1-page references + 7-page appendix. Dataset counts (Paper Table 1) and several COVID/H3 coverage cells match the locked e-abs protocol. The scientific story does not. Main-text Table 2 claims TreeSBM wins H3N2 RF / Branch-W1 / terminal edit with numbers that **do not exist** in any aggregate CSV; the same paper’s Appendix C.1 shows TreeSBM *worse* than Neutral CTMC on RF and W1. Table 4 reports COVID “key mutation recall” **0.219**, which is T8 `site_recall` (0.2189), not `mut_recovery` (0.080) — and C.2 correctly prints 0.080. H1N1 coverage in Table 4 is ε=3 while COVID/H3 are ε=2, and H1 min-dist is from the *old* morning job. Appendix C.3 relabels empty locked-root terciles as “short” and column-shifts the HIV-long row. Coverage@100 ≈ 0.775 is saturated on the five Brazil roots for almost every ablation; Tree-KL = 0.693 is ln2 everywhere. Fitness tilting Q0^F is written as part of the method, but the primary `covid_v5_mutrec` runs use **β=0**. Antibody results exist in the repo (and Neutral/CoSiNE beat TreeSBM on CDR) but are omitted from the PDF, while the intro still names antibody maturation as a motivating domain. A NeurIPS reviewer with the appendix open would treat Table 2 / Table 4 mut-recall as a reliability problem, not a presentation nit.

**Severity**

- **Critical**
  - Paper Table 2 TreeSBM RF/Quartet/W1/terminal-edit (0.9431 / 0.6336 / 4.45e-4 / 0.0067) vs repo TreeSBM **0.9841 / 0.7179 / 8.64e-4 / 0.0136** (`PAPER_TABLES.md`; `table_empirical_main.tex` RF **0.984 ± 0.002**). Caption “TreeSBM obtains the best RF” is false under repo numbers. PhylaFlow cells (0.9629 / 0.6210 / 1.08e-3 / 0.0083) have **no matching aggregate**.
  - Paper Table 4 COVID TreeSBM key-mut **0.219** vs repo **0.080** (`table5_viral_forecasting.csv` / T8 full `mut_recovery=0.0796`). 0.219 = `site_recall=0.2189`. Internally contradicted by Paper C.2 Mut **0.080**.
  - Paper C.3: COVID/H3 “short” cells are repo **medium**; HIV-long Cov/Mut are repo Mut/Clade (**column shift**). Empty locked buckets hidden. Conflicts with `table_c3_horizon.tex`.
- **Major**
  - H1N1 Table 4 Cov@100 **0.813 / 0.813 / 0.825** is `coverage_obs_e3` (job **7575016**); COVID/H3 use **e2**. H1 min-dist **4.04 / 3.88 / 3.43** is morning job **7539948**, not eabs **1.462 / 1.462 / 1.388**.
  - Table 2 vs C.1: same “TreeSBM on H3N2” cannot be both RF=0.9431 (best) and RF=0.983 (worse than Neutral 0.980).
  - Tree-KL 0.693 pasted as a ranking metric in C.1 / E.1 / E.2 / F.1 with no ln2 footnote; all methods tied.
  - Cov@100 flat **0.775** across Neutral / pLM±fit / all T8 ablations on locked Brazil roots; paper still frames coverage as a TreeSBM advantage (COVID 0.788 vs 0.775).
  - Q0^F / β>0 written as the method; primary viral ckpt and T8 full are **β=0**. Paper E.1 “TreeSBM default / fitness yes” is a **different** run (β=1, mut 0.038 vs Table 4 mut 0.080).
  - Paper E.1 caption: “TreeSBM uses each reference inside the same bridge-matching framework.” Repo T7 rows 1–6 are `--ablate-bridge` (pure R0).
  - Paper E.2 last two rows have **identical flags** (pLM+fit+branching) but AA|hit 0.086 vs 0.375; the real difference is **bridge**, which is not a column.
  - Antibody omitted from the PDF; repo Track C CDR: Neutral SHM **0.688** / CoSiNE **0.658** / pathogen TreeSBM **0.278** / OAS TreeSBM **0.150**. Intro still lists antibody maturation.
  - COVID pLM EVEscape Δ still **—** (job **7644883** PENDING). Paper correctly dashes it; do not fill.
- **Minor**
  - C.1 COVID Neutral pLM NLL paper **0.457** vs `table_c1_full_tree_gen.tex` **0.475** (copy of TreeSBM NLL). H3 AR/TreeSBM NLL swapped (0.463/0.462 vs tex 0.462/0.463).
  - Paper E.1 ProGen2 Mut **0.012** vs repo **0.1126**; ProGen2 Cov paper 0.775 vs **0.7875**.
  - Paper F.1 no-seq-branch uses old enrich (mut 0.045 / cons 0.824) not re-enrich paste (0.0375 / 0.818).
  - Figure numbering: body cites “Figure 4” then “middle panel of Figure 3”; PDF has Figures 1–3 + D.1 only.
  - EVEscape Δ direction: Table 4 **↑**, Appendix B.3 **lower**.
  - G.2 omits stop-threshold / seed / `branch_rate_scale=6.0` that `table_j2_j3_hyperparams.tex` has.
  - “North American H3N2” / “2023 and prior” vs collection-year 2024 test, train ≤2022 (`TABLE_METADATA.md`).
  - Grammar: “mutations occur are known antigenic sites”; “hyper-relevant.”

---

## 2. Methodology

### 2.1 Split protocols — what the paper says vs what was run

| Dataset | Paper wording | What the repo actually evaluated | Soundness |
|---|---|---|---|
| SARS-CoV-2 Spike | Table 1: **Geographic**. §3.4: “unseen region.” Table 3-style language (“future lineages”) is reused for COVID. | `prepare_covid_geo.py`: test=**Brazil** (40 trees, 2020-02-28…2025-09-09), val=Australia. **Not** `data/covid_temporal` (train≤2022 / test 2024–2025). Coverage roots: Brazil groups **4, 5, 7, 9, 10** (2021-05–2021-07). | Geographic holdout is a real test of spatial generalization, but it is **not** prospective forecasting. Brazil 2021 roots sit inside a pandemic that is already in the train countries. Do not call this “sequences originating before t.” |
| H3N2 HA | Table 1: **Temporal (pre-2024)**. Table 2: “held-out **2024 North American**.” Table 3: “training on sequences observed in **2023 and prior**.” | `prepare_h3n2_temporal.py`: train **≤2022** / val **2023** / test **2024 calendar year** (not NH Oct–Sep season; not “North American” as a geographic split). Coverage: 5 roots, all **group_002**, dates 2024-01-03…01-05 (`TABLE_METADATA.md`). | Closest thing in the paper to true forecasting. Two problems: (i) “2023 and prior” is wrong (train ends 2022); (ii) five roots from **one January 2024 group** is a narrow 2024 draw, not a season. |
| H1N1 HA | Table 1: **Geographic**. Table 4 dumps it next to COVID as the same “viral evolution forecasting” protocol. | `prepare_h1n1_geo.py` location-unit holdout (`data/h1n1`). Temporal 2025 set exists and was **not** used for `h1n1_v2_lit_hotspot` / Table 5. | OK if labeled geographic. Mixing H1 geo with H3 temporal under one “forecasting” heading overclaims time-travel. |
| HIV Env | Appendix A.1: Geographic, 147 trees / 43,363 seqs. C.1–C.3 unlabeled “HIV Env.” | Repo has **geo and temporal**. C.1/C.2 paste is **geo** (TreeSBM n_roots=8/14). Temporal RF/Q slightly better, W1 worse, n=6/14. | Choosing geo is fine; not naming it invites mixup (and C.3 already mixed columns). `table2_datasets.md` still says “HIV Env: skipped” for main Table 1 — HIV appears only in the appendix. |
| Antibody | Intro lists “antibody maturation.” No main or appendix results table in this PDF. | Three tracks (`antibody_tracks.md`): **A** OAS 1M donor-disjoint clones (`ab_oas_1m_v1`); **B** collapsed DASM `ab_dasm_v1`; **C** Rodriguez 82, **pathogen** `best.pt` (job **7062213** / **7517210**), forced observed topo+BL. Paper Track C Cov@100 uses **ε=5**, viral uses **ε=2**. | Omitting Ab results is safer than mixing tracks, but the intro then overclaims the domain. If Ab is added later: never label Thrifty as Neutral SHM; never put OAS and pathogen TreeSBM in one unlabeled “TreeSBM” row. |

**Clade holdout** (`data/covid_cladeholdout`) is documented in `SPLITS.md` and is **not** in this PDF. Do not imply clade generalization.

**Locked coverage roots.** COVID/H1/T7/T8/D.*/E.* share five roots. That is good for comparability and bad for variance: one hard COVID root (group 9, mean edit ≈5.6) sets the ε=2 plateau. Any claim that “coverage is robust across methods” is mostly “the same five trees.”

### 2.2 Coverage definition — are “coverage” conclusions valid?

Code (`COVERAGE_WHAT_MOVES_IT.md`): Cov@K,ε = fraction of **observed leaves** with min Hamming ≤ **ε absolute AA edits** to some generated leaf in a pool of K trees. Paper Table 3 caption (“ε = 2 mutations”) matches this. Legacy fractional `eps_frac=0.02` saturates at 1.0 on Spike L=1280 and is **not** the paper column — good that the main H3/COVID numbers use e-abs.

**Saturation on COVID Brazil.** NeutralBD, pLM±fit, AR, TreeSBM full, and every Table-8 ablation report `coverage_obs_e2 ≈ 0.775` at K=100 on groups 4/5/7/9/10. Easy leaves are within 0–2 AA of *some* gen leaf even for NeutralBD; group 9 is uncovered at ε=2 for everyone. K=10→100 barely moves e1/e2. **ε, virus, and genetic horizon move coverage; bridge / entropy / mask / pLM±fit do not.**

Therefore:

- Highlighting COVID Cov 0.788 vs 0.775 as a TreeSBM win is not a valid model comparison. 0.788 is within the same plateau (no_lit_mask is also 0.788).
- Ablation tables that show Cov@100 = 0.775 for every row cannot support “bridge matching improves coverage.”
- H3N2 e2 **does** move (Neutral 0.662 / TreeSBM 0.725). That is the only main-text coverage contrast with a real gap. Even there, Table 4-style ε=5 is flatter (0.800 vs 0.812).
- HIV geo e2 ≈ 0.188 for pLM and TreeSBM; temporal e2 = **0** (mean edit ~180). Coverage is a diversity/horizon metric, not a method ranking, once edit distances leave the ε ball.

Algorithm 6 in the PDF allows “edit distance, mutation-set distance, or pLM embedding distance.” The numbers in the tables are **AA Hamming**, not embeddings. Say so.

Paper “Median min dist” vs code `mean_min_edit` (mean, not median) — `table4_eps5.md`. Minor labeling, but it is a metric mismatch.

### 2.3 Tree-KL / Split-KL — misused as a ranking metric

Empirical Tree-KL/Split-KL are **NaN** in the CSVs (Dirac / unmatched clade sets). Values pasted as **0.693** are **sim_neutral Tree-JS ≈ ln 2 = 0.693147**, saturated for **every** method in jobs **7575017** / T8b **7585976–80**. Split-KL on sim_neutral *does* separate Neutral (~43) / TreeSBM (~57) / AR (~63) — and **Neutral wins**. Paper C.1 prints all three Tree-KL cells as 0.693 with ↓ arrows and no footnote.

Paper E.1 fills Tree-KL=0.693 for JTT / ESM-2 / ESM-C / ProGen2 even though repo D.1 marks per-R0 KL as **—** (not computed; would saturate). That is not “TreeSBM default happens to match ln2”; it is copying a vacuous constant into a comparison table.

**Reviewer attack:** “Your topology-distribution metric is identical for Neutral CTMC and TreeSBM. Why is it in every ablation table?” Honest move: drop Tree-KL from ranking tables; keep RF / Quartet / Branch W1 on empirical; optionally Split-KL on sim with an explicit “sim_neutral, not empirical” tag.

On empirical RF/W1, Neutral CTMC is as good or better than TreeSBM in C.1 (COVID RF 0.980 vs 0.982; H3 0.980 vs 0.983; HIV geo 0.972 vs 0.977; Neutral W1 strictly better). The method’s topology heads are **not** empirically justified by these tables. Contribution (1)–(2) “tree realism” is the claim most at odds with C.1.

### 2.4 Q0 vs Q0^F, β=0, ESM conservation prior

Paper eqs. (3)–(4): mutation prior Q0 = p_pLM(x′|x) (ESM-2 8M), then fitness tilt Q0^F ∝ Q0 exp(β F(x′)). The text presents tilting as part of the **reference process**, not an optional ablation.

**Implemented primary path** (`TABLE_METADATA.md`, T8 full, `covid_v5_mutrec`): `fitness_beta` unset → **β=0**. Table 4 TreeSBM mut **0.080** / cons **0.839** / aa|hit **0.375** is this run.

**β=1 run** is repo Table 7 “TreeSBM default (fitness β=1)”: mut **0.0378**, cons **0.941**, aa|hit **0.290**, Cov still 0.775. Fitness *hurts* mutation recovery (same pattern as ESM-2-650M fit: mut 0.001, cons 0.988). Paper E.1 still marks “TreeSBM default / Fitness weight **yes**” and quotes the β=1 mut 0.038, then Table 4 quotes the β=0 mut (or the site_recall mixup).

ESM-2 8M MLM is a **conservation** prior: it prefers plausible, often less-mutated sequences. That is why Cons is high for pLM/AR (0.985 / 0.993) and why TreeSBM’s mut-recovery gains come with **lower** Cons (0.839). The paper lists Cons with ↑ and does not discuss the tradeoff.

For antibodies, Q0 = ESM MLM is the wrong biological prior (SHM targets CDR/AID hotspots, not conserved framework). Repo `antibody_tracks.md` is explicit: v1 CDR=0.150 is a **Q0 mismatch**, not a missing viral mut-recovery flag; v2 SHM site-boost is queued and **not** in the PDF. If Ab stays out of the PDF, do not promise the domain. If Ab goes in, Q0^F with β>0 would likely **hurt** CDR further (anti-SHM).

### 2.5 Neutral SHM vs Thrifty; pathogen vs OAS TreeSBM

Not in this PDF (Ab table omitted) — still a landmine if a row is pasted from an older ICLR draft:

- **Neutral SHM** = JC69 independent-site NT CTMC on observed topo+BL. CDR **0.688**.
- **Thrifty** job **7517208** = additive context-SHM (`ThriftyHumV0.2-59`). CDR **0.746**. **Not** Neutral.
- Paper Track C **TreeSBM** = pathogen `checkpoints/best.pt` (flu, job **7062213**), forced topo. CDR **0.278**, Cov@e5 **0.084**.
- OAS `treesbm_ab_oas` (`ab_oas_1m_v1`, **7627919**) CDR **0.150**, Cov@e2 **0.024** (C.2 uses 0.027 @ 100 rollouts). Worse CDR than Neutral/CoSiNE; slightly better neighborhood coverage than Neutral.

Forced-topo Track C does **not** test tree generation (RF NaN). C.1 Ab free-topo (**7635753**) is a different experiment (Tree-KL again 0.693; TreeSBM W1 **worse** than Neutral). Do not mix.

### 2.6 C.3 locked roots / empty terciles

Protocol (`table_c3_horizon_protocol.md`): reuse the **same** Table 5 / Table 4 roots; assign each root to a train-H tercile; regenerate at that root’s own mean RTT H. Empty buckets (COVID short, H3 short, H1 med, HIV geo med, HIV temporal short) mean **none of the locked roots landed there**. Filling them would require new roots — forbidden.

**Honest presentation:** keep `--` and one sentence: “empty = no locked root in that tercile, not a missing job.”

**This PDF:** drops the `--` rows, labels COVID/H3 **medium** numbers as **short**, and shifts HIV-long (Cov 0.000 / Mut 0.050 / Clade 0.097 → paper Cov 0.050 / Mut 0.097). That looks like missing data *and* a spreadsheet error. `table_c3_horizon.tex` already has the correct `--` / medium / long layout (except COVID-long Mut 0.003, which is NeutralBD, not TreeSBM — even the tex is slightly wrong; TreeSBM long Mut is **0.000**).

H1N1 C.3, where TreeSBM actually looks good (short mut 0.212 vs Neutral 0.045; long Cov 0.500 vs 0.406), is **omitted**.

### 2.7 Ablations that do not move Cov@100 — narrative coherence

Paper F.1 (repo T8) is the right *experiment* (bridge / tree-context / BL head / entropy) and the wrong *lead metric*. Mut recall and aa|hit **do** move:

| Variant | Cov@e2 | Mut | AA\|hit | Cons | Source |
|---|---:|---:|---:|---:|---|
| w/o bridge | 0.775 | **0.008** | **0.086** | 0.876 | `table8_ablations.md` |
| full (β=0) | 0.775 | **0.080** | **0.375** | 0.839 | same |
| w/o entropy | 0.775 | 0.084 | **0.312** | 0.835 | entropy helps aa\|hit, not mut count |
| no_lit_mask | **0.788** | 0.069 | 0.341 | 0.839 | gen-time mask drop; v5 trained **without** PMC mask — not a train ablation |

Coherent narrative: **bridge matching is about mutation identity (mut / aa|hit), not set-cover at ε=2.** Incoherent narrative: “full TreeSBM improves coverage” when every row is 0.775.

no_bl_eval ≈ full on mut/aa|hit (0.080 / 0.375): the BL head is not doing work at gen time on this protocol. no_tree_ctx mut 0.084 > full 0.080: tree context is not a mut-recovery win here. Say that; do not rank by Tree-KL=0.693.

**Sampling vs training ablations.** T8 flags (`--ablate-bridge`, `--ablate-tree-context`, …) are **eval-time** on `covid_v5_mutrec`. Repo E.1 terminal-only / no-Doob / no-terminal-consistency are **retrains**. Paper F.1 is the eval-flag table; Paper E.2 “pLM + branching” vs “full bridge” confuses eval-flag factorial with a trained component ablation. Caption must say “frozen ckpt, generation-time ablate” vs “new train.” λ_br=0 retrain **7539981** TIMED OUT — do not use.

### 2.8 Fairness of baselines

| Baseline | What it is in code | Fairness issues |
|---|---|---|
| Neutral CTMC + BD | Neutral substitution + birth-death; often forced or simple BD topology | Strongest topology baseline in C.1. Paper Table 2 still ranks it below a TreeSBM row that is not in the CSV. Underclaiming Neutral on trees is the opposite of the usual baseline-starving problem. |
| pLM mutation prior | ESM-2 8M on a tree (observed or simple BD) | Fair sequence-only control. COVID EVEscape for this row is **NaN** (ESM L=566 OOB on L=1280, job **7539949**); Table 4 dashes it — correct. Do not compare EVEscape only on AR/TreeSBM and imply pLM was scored. |
| Autoregressive tree-edit | **ARTreeFormer-adapted pool**, N=16 prune + **JTT**, not an Ab-trained or virus-specific AR | Adapted-pool, not native forward gen (`BLOCKERS.md`). Ab AR CDR=**0.000** is a real JTT outcome, not a stub. Calling it “autoregressive tree-edit” without “adapted pool” overstates the competitor. |
| JTT + BD / EmpiricalBD | “not literal JTT CTMC — caption carefully” (`PAPER_TABLES.md`) | Paper Table 2 name is acceptable if caption says empirical substitution + BD. |
| PhylaFlow | Wired; N=16 empirical RF **0.987**, Quartet **0.672** (`table_empirical_N16.csv`) | Paper’s 0.9629 / 0.6210 / TE 0.0083 is **not** that aggregate. Native vs adapted must be labeled (`EXTERNAL.md`). |
| CoSiNE / Thrifty / DASM | Strong Ab CTMC baselines | Absent from related work and from the PDF. If Ab is in scope, CoSiNE is the obvious reviewer question. |
| TreeSBM vs pLM on coverage | Same 5 roots, same ε | Fair protocol, weak discriminator. Mut/aa|hit is the fair fight; TreeSBM does win those on COVID (0.080 / 0.375 vs pLM 0.048 / 0.017) **if you use mut_recovery**. |

Enrichment (n_trees=20, mrs=0.5, n_steps=100) vs coverage (K=100, n_steps=50, seed 0 vs enrich seed 42) are **different generators**. Mixing coverage-CSV `mut_recovery` (Table 3 TreeSBM **0.156**) with enrich `mut_recovery` (C.2 **0.100**) is a protocol collision, not a rounding issue.

---

## 3. Coherence

### 3.1 Abstract / intro / contributions vs tables

**Abstract:** “realistic evolutionary tree generation, future-lineage forecasting, and probabilistic exploration.”  
**Contribution (2):** “advantages in tree realism, future-sequence coverage, and mutation recovery.”  
**Contribution (3):** “generated mutations can recover known escape-associated variation and pandemic-grade VOCs.”

| Claim | Evidence in this PDF + repo | Verdict |
|---|---|---|
| Tree realism | C.1: Neutral ≤ TreeSBM on RF; Neutral **better** W1. Table 2 TreeSBM win uses unaudited numbers that contradict C.1. | **Not supported** by authoritative artifacts. |
| Coverage advantage | H3 e2: 0.725 vs Neutral 0.662 — real, modest, 5 roots from one group. COVID e2: 0.788 vs 0.775 — saturated. H1 e2: **tied 0.762** (paper instead reports e3). HIV geo: tied 0.188. | **Overclaimed** as a general forecasting win. |
| Mutation recovery | COVID enrich mut 0.080 vs pLM 0.048 / AR 0.006; aa|hit 0.375 vs 0.017 / 0.094. H1 mut 0.241 vs 0.071. H3 enrich 0.100 vs coverage-CSV 0.156 vs pLM 0.004. | **Supported** if you quote **mut_recovery / aa|hit**, not 0.219, and admit Cons is worse. |
| Antigenic / VOC | COVID antigenic TreeSBM **0.005** vs AR **0.017** (TreeSBM **loses**). EVEscape Δ TreeSBM **−0.048** (generated less escape-like than GT). Figure 2/3 is qualitative Gamma recovery on one São Paulo 2020 root. | **Under-specified.** Qualitative VOC walk is interesting; Table 4 antigenic/EVEscape do not show TreeSBM beating AR on antigenic, and Δ is negative. B.3 says Δ lower-is-better; Table 4 marks ↑. |
| Antibody / protein families / directed evolution | Named in intro. Table 1 has three viral datasets only. No PFAM jobs (`APPENDIX_COMPLETION_PLAN.md`). | **Domain overclaim.** |

### 3.2 Antibody story vs viral story

Viral: TreeSBM can beat pLM on **which** amino acid is recovered given a hit, at the cost of conservation.  
Antibody (repo, not in PDF): Neutral JC69 and CoSiNE recover **CDR mutations** far better than TreeSBM (0.688 / 0.658 vs 0.278 pathogen / 0.150 OAS). That is expected if Q0 is ESM conservation. Putting both stories in one paper without the Q0 discussion looks like two different methods.

If the camera-ready stays viral-only, delete antibody/protein-family sentences from the intro. If Ab is added, the honest sentence is: “On affinity-maturation families, a JC69/CoSiNE SHM process matches CDR better than ESM-grounded TreeSBM; TreeSBM is not yet an SHM model.”

### 3.3 “Best model” in one table, not another

- Table 2: TreeSBM best RF/W1/TE (unaudited). C.1: Neutral best RF/W1. **Contradiction.**
- Table 3: TreeSBM best Cov@100 and min-dist; AR best clade recall (**0.300** vs TreeSBM **0.100**). Caption does not mention AR winning clades.
- Table 4: TreeSBM best mut/aa|hit/min-dist/EVEscape-closer-to-zero; AR best COVID antigenic and Cons. “Substantially higher key-mutation recall” is true only for mut_recovery, and the printed 0.219 is the wrong statistic.
- E.1: ProGen2 mut **0.113** (repo) would **beat** TreeSBM β=1 mut 0.038 and is competitive with β=0 mut 0.080, with Cons 0.986. Paper prints ProGen2 mut **0.012**, which hides that. Fitness-on ESM-2-650M has the best Cons and worst mut — the interesting result — but the caption still frames “TreeSBM default” as the reference-process winner.

### 3.4 Appendix vs main (ckpt, ε, N, virus)

| Item | Main | Appendix / repo | Issue |
|---|---|---|---|
| H3 TreeSBM ckpt | Table 2 topology historically `h3n2_v2`; Table 3 coverage `h3n2_v3_lit_hotspot` | Never named in the PDF | Two checkpoints, one name “TreeSBM.” |
| COVID TreeSBM | Table 4 = `covid_v5_mutrec` β=0, job **7575015** / T8 full | E.1 “default” = β=1 (`table7_treesbm_fit`) | Same name, different β. |
| ε | Table 3: ε=2. Table 4: unspecified. | COVID/H3/C.2/C.3: e2. H1 Table 4: **e3**. Ab (if added): e5. | Unlabeled ε change. |
| N | Unstated in main tables | N=16 everywhere for coverage | State N=16, K=100, 5 roots. |
| HIV | Absent in Table 1 | C.1–C.3 “HIV Env” = **geo** | Name the split. |
| Mut recall H3 | Table 3: **0.156** (coverage CSV @K=100) | C.2: **0.100** (enrich mrs=0.5) | Two protocols, one name. |

Best-ckpt identity the rest of the repo treats as locked: COVID **`covid_v5_mutrec`** (not v7; morning coverage **7539947** used v7 + degenerate group_001 — do not revive). H3 coverage **`h3n2_v3_lit_hotspot`**. H1 **`h1n1_v2_lit_hotspot`**. Paper never writes these names.

### 3.5 Notation vs implementation

- Eq. (7)–(8): log R_θ = log R0 + c_θ. Matches code (`RateHeads`, Doob-style additive control).
- Algorithm 1 matches the GSBM-style *matching* loop (`GSBM_BRIDGE_REVIEW.md`): single-stage, fixed empirical coupling, **no** CondSOC / IPF. Calling it a “tree-valued Schrödinger bridge” is a modeling slogan; the trained object is **bridge matching to a hand-specified conditional path** (`SampleBridgeState` time-cut + Bernoulli edits). A theory-aware reviewer will ask for the gap vs De Bortoli / Shi / Liu GSBM.
- Algorithm 4 samples a stop event; J.3 / code: **stop head is not used at inference** (fixed `n_steps`). Paper G.2 omits that; `table_j2_j3_hyperparams.tex` says “Stop threshold none.”
- Q0^F in §2.2 vs β=0 in the winning Table 4 row: the notation is implemented as a **switch**, not the default.
- “pLM site score” F(x′) is not defined at the level needed to reproduce (ESM PLL? mean token NLL? EVEscape?).

### 3.6 Related work / CoSiNE positioning

Cited: ProtTrans, ESM, ProGen2, ProteinGym, FastTree, IQ-TREE, RAxML-NG, BEAST, Hie 2021, EVEscape (Thadani 2023), DDPM, flow matching, D3PM, Léonard, DSBM, DSBM-matching.

**Missing relative to the experiments you actually run:** CoSiNE, Thrifty/netam, DASM, ARTreeFormer (only a table row), PhylaFlow (table row, no citation block), PhyloVAE, Neutral BD as a literature object, Nextstrain/TreeTime (methods A.1 names them). Schrödinger-bridge matching citations are present; GSBM (Liu et al.) is the closer matching ancestor and is uncited (`GSBM_BRIDGE_REVIEW.md`).

CoSiNE is the standard neural CTMC for affinity maturation. Silence looks like unawareness, not a scoping choice, because the intro already invoked antibodies.

---

## 4. Result consistency audit

For every main + appendix table in the PDF. Paper value quoted from pypdf extract; repo from the cited file.

### Paper Table 1 — Datasets (p.4)

| Cell | Paper | Repo (`table2_datasets.md`) | Status |
|---|---|---|---|
| COVID trees / seqs / split | 404 / 104,260 / Geographic | 404 / 104,260 leaves / Brazil geo | **Match** |
| H3N2 | 295 / 118,000 / Temporal (pre-2024) | 295 / 118,000 / train≤2022, test 2024 | **Match counts**; “pre-2024” OK; not “North American” |
| H1N1 | 301 / 55,755 / Geographic | 301 / 55,755 / geo location-unit | **Match** |
| HIV / Ab / PFAM | absent | HIV appendix-only; Ab 82 / 377 Rod.82; PFAM none | OK for this table; intro still names extra domains |

### Paper Table 2 — H3N2 topology (p.4)

Baseline rows match `PAPER_TABLES.md` TEMP empirical (n_roots=97): Neutral 0.9824 / 0.7110 / 7.22e-4 / 0.0148; JTT 0.9839 / 0.7119 / 7.22e-4 / 0.0151; pLM 0.9823 / 0.7128 / 7.22e-4 / 0.0151; ART 0.9891 / 0.7126 / 8.61e-4 / 0.0151; PhyloVAE 0.9900 / 0.7160 / 8.12e-4 / 0.0151.

| Method | Paper RF / Q / W1 / TE | Repo | Status |
|---|---|---|---|
| TreeSBM | **0.9431 / 0.6336 / 4.45e-4 / 0.0067** | **0.9841 / 0.7179 / 8.64e-4 / 0.0136** (`PAPER_TABLES.md`); tex **0.984 ± 0.002 / 0.718 ± 0.006 / 0.001 / 0.014**; C.1 N=16 **0.983 / 0.721 / 9.25e-4** | **Critical — not in any aggregate CSV** |
| PhylaFlow | **0.9629 / 0.6210 / 1.08e-3 / 0.0083** | N=16 emp **0.987 / 0.672 / 9.23e-5 / 0.0055** (`table_empirical_N16.csv`) | **Critical — unmatched** |

Internal: C.1 TreeSBM RF **0.983** vs Table 2 **0.9431**. Caption “best RF / W1 / terminal edit” is false under repo numbers (Neutral RF 0.9824 < TreeSBM 0.9841).

### Paper Table 3 — H3N2 future lineage (p.4)

ε=2 caption is correct. Cov@100 matches e2 (`table4_coverage_vs_eps.csv` / `TABLE_METADATA.md` companion): Neutral 0.662, pLM 0.675, AR 0.675, TreeSBM **0.725**. Text “72.5% … within two mutations” matches.

| Cell | Paper | Repo | Status |
|---|---|---|---|
| Neutral Cov@10 / @100 | 0.662 / 0.662 | e2 K=10 **0.6625** / K=100 **0.6625** | Match (round) |
| pLM Cov@10 / @100 | 0.662 / 0.675 | e2 0.6625 / **0.675** | Match |
| AR Cov@10 / @100 | **0.622** / 0.675 | e2 **0.6625** / 0.675 | **AR@10 mismatch** (likely typo 0.662→0.622) |
| TreeSBM Cov@10 / @100 | 0.675 / 0.725 | e2 **0.675** / **0.725** | Match |
| Mut recall N / pLM / AR / TS | 0.020 / 0.004 / 0.079 / **0.156** | coverage CSV K=100 (`table4_eps5.md`) **0.156**; C.2 enrich **0.100** | Table 3 vs C.2 **internal protocol clash** |
| Clade AR / TS | 0.300 / 0.100 | `table4_eps5.md` 0.300 / 0.100 | Match (coverage path) |
| Min dist | 2.65 / 2.65 / 2.40 / 2.24 | same | Match (`mean_min_edit`, not median) |

ε=5 companion (job **7459412**) is 0.800/0.812 — do not silently switch.

### Paper Table 4 — Viral forecasting (p.5)

COVID coverage/aa|hit/min-dist/EVEscape/cons (except mut) match `table5_viral_forecasting.csv` job **7575015**, ckpt `covid_v5_mutrec`. pLM EVEscape dashed — matches **7644883** PENDING. Cons paper 0.840 vs 0.8389 — OK.

| Cell | Paper | Repo | Status |
|---|---|---|---|
| COVID TSBM key mut | **0.219** | **0.0796** mut_recovery; **0.2189** site_recall (`table8_ablations.md` full) | **Critical — wrong statistic** |
| COVID TSBM antigenic | 0.005 | 0.0051 | Match; **worse than AR 0.017** |
| H1 Cov@100 pLM / AR / TS | **0.813 / 0.813 / 0.825** | e2 all **0.7625**; **e3** 0.8125 / 0.8125 / **0.825** (`coverage_curves_h1n1_N16_table5_eabs.csv` K=100) | **Major — ε=3 vs COVID ε=2** |
| H1 min dist | **4.04 / 3.88 / 3.43** | eabs **1.462 / 1.462 / 1.388**; old **7539948** min_edit **4.04 / 3.88 / 3.43** | **Major — mixed jobs** |
| H1 mut / aa|hit / antigenic / cons | 0.071 / 0.006 / 0.241; … 0.618; 0.163; 0.945/0.998/0.952 | `table5_viral_forecasting.csv` | Match (enrich) |
| H1 AR EVEscape | −0.297 | −0.2975 | Match (round) |

### Paper A.1 (p.8)

HIV 147 / 43,363 Geographic — not in `table2_datasets.md` main paste (HIV skipped there). A.1 construction (MAFFT, FastTree, Augur, TreeTime) matches the pipeline description. Caption still says splits may hold out “protein families, clades, or full roots” — those evals are **not** in this PDF.

### Paper B.1–B.3 (p.8–9)

Glossary is directionally right. Conflicts: B.3 EVEscape **lower** vs Table 4 **↑**; “Median min distance” vs mean; Coverage@K,ε does not specify absolute vs % (Table 3 does). Tree-KL listed as if informative.

### Paper C.1 (p.9) vs `table_c1_full_tree_gen.tex`

Topo RF/Q/W1/Split-KL match the tex for COVID / H3 / HIV **geo**. Tree-KL 0.693 all rows — saturated, no †.

| Cell | Paper | Tex / plan | Status |
|---|---|---|---|
| COVID Neutral NLL | **0.457** | tex **0.475**; md had — then job **7635749** | **Copy-paste of TreeSBM 0.457** |
| H3 AR / TS NLL | **0.463 / 0.462** | tex **0.462 / 0.463** (`0.4627` job **7635752**) | **Swapped** |
| HIV | geo numbers, unlabeled | geo n=8/14 | OK numbers; name geo |
| Antibody rows | **missing** | tex filled (NLL 0.513 / 0.466 / 0.451) | Appendix incomplete vs tex |
| Protein families | missing | — no jobs | OK to omit |

### Paper C.2 (p.9) vs `table_c2_future_lineage.tex`

Viral six rows **match** the tex, including COVID TSBM Cov@500 **–** (job **7635911** RUNNING, CSV 0 bytes). Antibody two rows in tex (**pLM 0.018 / TSBM 0.027**) **dropped from PDF**. HIV unlabeled geo. Mut COVID TSBM 0.080 **contradicts Table 4 0.219**. Mut H3 TSBM 0.100 **contradicts Table 3 0.156**.

### Paper C.3 (p.10) vs `table_c3_horizon.md` / `.tex`

| Paper row | Paper numbers | Repo | Status |
|---|---|---|---|
| COVID “short” 0.938 / 0.00 / 0.69 | medium TreeSBM (all methods tied) | COVID **short = empty** | **Critical mislabel** |
| COVID long 0.667 / **0.003** / 1.92 | NeutralBD mut **0.003**; TreeSBM mut **0.000**, min 1.92 | Tex also has 0.003 | Mut is Neutral, not TreeSBM |
| H3 “short” 0.859 / 0.099 / 1.14 | **medium** TreeSBM | H3 short empty | **Critical mislabel** |
| H3 long 0.00 / 0.131 / 7.19 | TreeSBM long | Match | OK |
| HIV short 0.313 / 0.001 / 88.6 | geo TreeSBM short | Match | OK |
| HIV long **0.050 / 0.097 / 24.7** | repo Cov **0.000**, Mut **0.050**, Clade **0.097**, min 24.7 | **Column shift** | **Critical** |
| H1N1 | absent | filled; TreeSBM strongest mut/cov among C.3 | Omitted evidence |

Empty buckets explained in protocol; PDF makes them look like short-horizon data.

### Paper Fig D.1 (p.10)

Qualitative Brazil Gamma 2020 / Rio Omicron 2023. Repo H.1 documents COVID g40 + H3 screens; HIV/Ab galleries **—**. Fine as illustration; not a result. Body “Figure 4” / “Figure 3” numbering is broken relative to Fig 2 middle (Gamma counts) vs Fig 3 (trajectory).

### Paper E.1 (p.10) ≈ repo D.1 / `table7_refproc.md`

Caption “each reference inside the same **bridge-matching** framework” is **false** for rows 1–6 (`ablate_bridge=true`). Tree-KL 0.693 filled for R0 backends (repo —).

| Row | Paper Mut / Cov / NLL | Repo Mut / Cov e2 / NLL | Status |
|---|---|---|---|
| JTT | 0.023 / 0.775 / 0.471 | 0.0223 / 0.775 / 0.4705 | OK |
| ESM2 nofit | 0.018 / 0.775 / **0.445** | 0.0180 / 0.775 / **0.4415** | NLL off |
| ESM2 fit | 0.001 / 0.775 / 0.441 | 0.0011 / 0.775 / 0.4405 | OK |
| ESM-C ±fit | 0.035 / 0.036; 0.461 / 0.427; Cov 0.775 | 0.0347 / 0.0356; 0.4605 / 0.4267; Cov **0.775** (jobs **7628260/61** COMPLETED) | OK; PASTE_TABLES still said Cov — — **stale paste vs table7** |
| ProGen2 | Mut **0.012**, Cov **0.775**, NLL 0.454 | Mut **0.1126**, Cov **0.7875**, NLL 0.4546 (**7627918/7628262**) | **Major** |
| TreeSBM default fit yes | 0.038 / 0.775 / 0.446 | β=1: 0.0378 / 0.775 / 0.4464 | Matches **β=1**, not Table 4 β=0 |
| Evo2 | absent | — backend not wired | OK |

### Paper E.2 (p.10) ≈ repo D.2

AA|hit 0.108 / 0.098 / 0.375 / 0.086 / 0.375 matches D.2 map. Cov all 0.775 — true and uninformative. Tree-KL 0.693 on ESM2 rows is repo —. Last two rows **same checkbox flags**, different AA|hit: real contrast is T8 full vs T8 **no_bridge**, not “fitness+branching vs full.” D.2 notes “pLM + branching” **is** T8 full (β=0). Paper implies a factorial that was not run.

### Paper F.1 (p.11) ≈ repo Table 8 / E.2 architecture

| Variant | Paper Mut / Cons / AA\|hit / NLL / Tree-KL | Repo paste | Status |
|---|---|---|---|
| w/o bridge | 0.008 / 0.876 / 0.087 / 0.447 / 0.693 | 0.0082 / 0.8758 / 0.0858 / 0.4465 / 0.6931 | OK |
| w/o seq-dep branch | **0.045 / 0.824 / 0.247** / 0.454 / **0.693** | re-enrich **0.0375 / 0.8181 / 0.2418**; Tree-KL **—** | Old JSON; KL invented |
| w/o BL / tree-ctx / entropy / full | ~paste | match (round) | OK; tree-ctx KL paper 0.693 vs repo — |
| Cov | all 0.775 | all 0.775 (no_lit_mask 0.788, not in F.1) | Honest numbers, dishonest ranking if used as a win |
| Fitness row | omitted | not a T8 row (β=0 already) | Correct omission |

Repo E.1 terminal-only / w/o terminal (Mut 0.089 / 0.062) and E.2 mut/stop-head trains are **not** in the PDF (no_doob still RUNNING).

### Paper G.1 / G.2 (p.11) vs `table_j2_j3_hyperparams.tex`

Layers 4, d=128, 8 heads, ESM-2 8M, head depth 2, dropout 0.1, AdamW, 1e-4, batch 1, epochs **14** — **match** `covid_v5_mutrec` job **7350224**. Horizon 50, max nodes 400, binary branching, T=1.0, K=100 — match Table 5 coverage protocol. Missing vs tex: stop none, seed 0, branch scale 6.0. Enrichment uses n_steps **100**, mrs **0.5** — not in G.2; say coverage vs enrich.

### Paper H Algorithms (p.12–14)

Useful. Gaps vs code: stop head unused at sample time; SampleBridgeState is time-cut + Bernoulli, not a learned Schrödinger potential; coverage alg allows embedding distance unused in tables. Algorithm numbering skips a “5” in the extract (4 then 6) — check the compiled PDF.

### Cross-cutting: ckpt / ε / saturation (requested checks)

- **Best COVID ckpt `covid_v5_mutrec`:** consistent for Table 4 / C.2 / T8 / F.1 / G.1. Broken by E.1 β=1 “default” and by any leftover v7 coverage.
- **ε:** H3 Table 3 and COVID Table 4 / C.2 = e2. H1 Table 4 = e3. Ab (repo) = e5. Fractional ε unused in the printed COVID/H3 cells — good.
- **Tree-KL ≈ ln2:** presented as discriminative in C.1/E/F — **not honest**.
- **Cov@100 = 0.775:** paper is **not** honest in F.1/E.2 that this is a plateau; slightly honest in leaving ProGen2/Evo gaps; dishonest in treating 0.788 as a method win.
- **C.3 empty buckets:** **not** explained; relabeled.
- **HIV geo vs temporal:** C.1–C.3 are geo; temporal exists and is omitted. C.3 long column-shift is worse than a silent geo choice.
- **COVID pLM EVEscape:** correctly **—**; job **7644883** still PENDING. Do not invent.

---

## 5. Scientific / writing feedback (reviewer-style)

### What a NeurIPS reviewer would attack

1. **Reliability of Table 2 / Table 4 mut 0.219 / C.3.** Appendix disagrees with the main text on the same datasets. This is a “can I trust any number?” review, not a taste comment.
2. **Is this a Schrödinger bridge?** Matching to a fixed tree-valued interpolation is closer to discrete bridge matching than to a computed Schrödinger potential. Related-work gap (GSBM) plus Algorithm 1 will be probed.
3. **Forecasting vs geography.** COVID and H1N1 are spatial holdouts; only H3N2 is temporal, and that test is five roots from one January group. “Evolution forecasting” as the title-level claim is ahead of the design.
4. **Topology.** If Neutral BD matches or beats TreeSBM on RF and W1, why the graph transformer and branching heads? Either show a split where topology heads matter (free-topo Ab C.1 does **not** help TreeSBM W1) or recast the paper as **sequence-on-a-tree / mutation-bridge**.
5. **Coverage as a KPI** with ε=2, K=100, n=5, N=16. Reviewers who plot the ε-sweep (`COVERAGE_WHAT_MOVES_IT.md`) will see the plateau.
6. **VOC / EVEscape.** Negative Δ and antigenic 0.005 vs AR 0.017 vs a figure that “recovers Gamma.” Pick one story: rare-path qualitative recovery, or quantitative escape enrichment — not both without reconciliation.
7. **Missing CoSiNE / CTMC literature** given intro antibodies and Neutral CTMC baselines.
8. **No error bars / n.** n_roots=5 for coverage; HIV TreeSBM n=8; Table 2 older n=97 vs C.1 N=16 means. NeurIPS will ask for uncertainty.

### Missing limitations (beyond p.6)

Page 6 mentions reconstructed trees, sampling, no recombination/indels. Add:

- Evaluation n is tiny for coverage (5 locked roots); one hard root dominates COVID ε=2.
- Tree-KL is non-informative on the empirical track.
- Primary model uses **β=0**; fitness tilt reduces mut recovery.
- Stop head trained, not used.
- AR / PhylaFlow are adapted-pool, not native generators.
- Antibody Q0 is ESM conservation; SHM not modeled (if Ab stays in the intro).
- Coverage ≠ antigenic forecasting; PMC hotspot fraction is small for all methods (~0.005).
- Temporal leakage risk: COVID Brazil test dates overlap the global pandemic years in train countries (geographic ≠ future).

### Overclaim / underclaim

**Overclaim:** tree realism; coverage as the headline metric; VOC recovery as a tabled result; Q0^F as the operating reference; “2023 and prior”; “North American”; antibody/protein-family scope.

**Underclaim:** the actual mut / aa|hit gap vs pLM on COVID and H1 (real, if you print 0.080 not 0.219); H1 C.3 horizon where TreeSBM mut/cov move; Cons–mut tradeoff (you have the numbers); Neutral CTMC as a strong topology baseline.

### Suggested rewrites (short)

**Abstract last sentence.** Replace forecasting boilerplate with: “On held-out viral trees, TreeSBM improves amino-acid identity at recovered mutation sites relative to an ESM mutation prior, with little change in ε=2 leaf coverage on a saturated SARS-CoV-2 split.”

**Contribution (2).** “On influenza HA and SARS-CoV-2 Spike, TreeSBM recovers more true amino-acid substitutions than pLM and AR baselines (aa|hit, mut recall). Topology distances remain comparable to a neutral birth-death process.”

**Table 2 caption.** Delete “best RF.” Either paste `table_empirical_main.tex` TreeSBM **0.984 / 0.718 / ~8.6e-4 / 0.014** or drop Table 2 until the 0.9431 artifact is found.

**Table 4.** Key mutation recall = `mut_recovery` **0.080**. State ε=2 absolute, N=16, 5 Brazil roots, ckpt `covid_v5_mutrec`, β=0. H1: same ε as COVID, min-dist from **7575016**. Footnote: Cov@ε=2 is saturated for this root lock.

**§2.2.** “We optionally tilt Q0 by a pLM fitness score (β>0). Unless noted, experiments use β=0; §E.1 reports β=1.”

**C.3.** Use `table_c3_horizon.tex` (keep `--`). One sentence on locked roots. Fix HIV-long. Use TreeSBM mut 0.000 for COVID long, or show all methods.

**F.1 caption.** “Coverage@100 (ε=2) is 0.775 for all rows on five locked roots and is not used for ranking. Bridge matching changes mut recall (0.008 → 0.080) and aa|hit (0.086 → 0.375).”

---

## 6. Appendix completeness vs `APPENDIX_COMPLETION_PLAN.md`

NeurIPS appendix is a **subset** of the ICLR plan, with **renumbered** sections. Remaining `--` / gaps:

| Plan ID | NeurIPS PDF | Status vs plan |
|---|---|---|
| A.1–A.2, B.* | A.1 + B.1–B.3 | A.2 split table not separate (folded into Table 1 / A.1). |
| C.1 | C.1 | Viral filled; **Ab tex rows dropped**; PFAM —; Tree-KL unfootnoted; Neutral COVID NLL wrong. |
| C.2 | C.2 | Viral match; **Ab rows dropped**; COVID TSBM Cov@500 still **—** (**7635911**). |
| C.3 | C.3 | Filled but **mislabeled / shifted**; H1 omitted; Tree-KL — (OK); empty buckets not shown as `--`. |
| D.1 pLM ref | Paper **E.1** | Mostly filled; ProGen2 mut wrong; Evo2 —; Tree-KL invented for R0. |
| D.2 components | Paper **E.2** | Mapped, not a real factorial; flags mislabeled. |
| E.1 bridge Doob/terminal | **absent** | Terminal-only / w/o terminal **filled in repo**; no_doob **RUNNING**. |
| E.2 mut/stop heads | **absent** | trains **7642530/31** queued; metrics — |
| E.3 sampling | **absent** | temp/maxl **filled** (flat 0.775); K=500 — |
| F.1–F.2 scaling/transfer | **absent** | Defer (no jobs) — OK |
| G.1 ranking / Fig G.1 | **absent** | Defer — OK |
| H.1–H.2 qualitative | Fig D.1 only | Partial |
| I.1 compute | **absent** | Repo filled: 11M, 8.0 h, 71 min/root; GPU mem — |
| J.1 probes | **absent** | FINAL in `PAPER_TABLES.md` (H3N2+COVID) |
| J.2–J.3 | G.1–G.2 | Present; slightly trimmed |
| K algorithms | H | Present |
| Main Ab table | **absent** | Repo Track C filled; Neutral ≠ Thrifty |

**Do not** stuff F/G/Evo2/C.3 Tree-KL before submit. **Do** fix C.3/C.1/E.1 numbers against existing tex, restore or drop Ab consistently, and add one coverage-saturation sentence.

---

## 7. Priority fix list (PDF before submit)

**Do first (wrong or contradictory numbers)**

1. Replace Table 2 TreeSBM **and** PhylaFlow with `table_empirical_main.tex` / C.1 H3 numbers, or remove PhylaFlow until an artifact exists. Rewrite the “best RF” sentence.
2. Table 4 COVID key-mut **0.219 → 0.080** (`mut_recovery`). Confirm you are not printing `site_recall`.
3. Rebuild C.3 from `table_c3_horizon.tex`: `--` for empty terciles; medium not short; HIV long Cov **0.000**, Mut **0.050**, Clade **0.097**. COVID long Mut TreeSBM **0.000** (not Neutral 0.003).
4. H1 Table 4: one ε (e2: Cov **0.762** all methods) and min-dist **1.46 / 1.46 / 1.39** from **7575016**. Do not mix **7539948**.
5. C.1 Neutral COVID NLL **0.475**; unswap H3 AR/TS NLL; footnote Tree-KL = ln2, not a rank.
6. E.1 ProGen2 mut **0.113**, Cov **0.788**; do not put Tree-KL on R0 rows; caption **ablate-bridge** for rows 1–6. Name β=1 vs Table 4 β=0.

**Do next (honesty / reviewer)**

7. One footnote on Cov@ε=2 Brazil plateau (`COVERAGE_WHAT_MOVES_IT.md`). Rank F.1 by mut / aa|hit.
8. State β=0 for `covid_v5_mutrec`; define Q0^F as optional.
9. Name ckpts: `covid_v5_mutrec`, `h3n2_v3_lit_hotspot` (coverage), `h3n2_v2` if Table 2 stays v2.
10. Fix Figure 3/4 numbering; EVEscape Δ arrow vs B.3; “median” vs mean min dist.
11. Intro: drop antibody/protein-family/directed-evolution or add a scoped sentence + CoSiNE citation. Do not paste Thrifty as Neutral SHM.
12. Table 3 vs C.2: pick enrich **or** coverage-CSV mut recall for H3 and use it in both places (0.100 vs 0.156).
13. AR Cov@10 0.622 → **0.662**.

**Optional if space**

14. I.1 compute (8.0 h, 11M) — reviewers ask.
15. H1 C.3 (actual horizon signal).
16. Restore C.1/C.2 antibody rows **only** with Track labels (pathogen vs OAS; ε=5 vs e2).
17. Leave **—**: Cov@500 COVID TSBM (**7635911**), COVID pLM EVEscape (**7644883**), Evo2, F/G, no_doob, E.2 mut/stop heads.

**Do not** invent any of those pending cells to “complete” the appendix.
