# Filling paper Tables 2 / 3 / 5 (skip Table 4)

**Canonical fill-in + remaining-work plan (all main + appendix tables, antibody §5.5/Table 6, TEMP vs FINAL):**  
[`checkpoints/ICLR_TABLE_FILL_PLAN.md`](../checkpoints/ICLR_TABLE_FILL_PLAN.md)  
(aligned to PDF `TreeSBM__ICLR_ (6).pdf`).

**Do not invent numbers.** Every cell must come from a JSON/CSV artifact produced by the commands below. Empty cells stay empty until the corresponding job finishes.

Splits: temporal forecasting + COVID clade holdout are documented in [`SPLITS.md`](SPLITS.md). Geographic dirs (`data/h3n2`, `data/h1n1`, `data/covid`) remain the default Wave-1 eval paths unless a caption says temporal/clade.

**TEMP vs FINAL:** cells below marked **TEMP** are draft placeholders only (mut_recovery still ≪ paper story; caption protocol must match artifact before FINAL). **H3N2 Coverage@K curves filled** (job **7439926** COMPLETED → `benchmarks/results/coverage_curves_h3n2_N16.csv` + `plots/coverage_curves_h3n2_N16.{png,pdf}`). COVID Track A Coverage@K still separate (7431995/96 TIMED OUT; see `results/wave5_pull/track_a/`).

---

## Wave-5 train status (2026-08-07) — see [`WAVE5_EVAL_RESULTS.md`](../checkpoints/WAVE5_EVAL_RESULTS.md)

| Job | Ckpt | State | Best ep / val | Enrichment |
|----:|------|-------|---------------|------------|
| **7434287** | `h3n2_v3_lit_hotspot` | **COMPLETED** | ep57 / 5.1519 | **RUNNING** 7436591 (mrs 0.3/0.5 + flu lit) |
| **7429545** | `covid_v5_deepmut` | **TIMEOUT** (24h; best kept) | ep10 / 2.8521 | **RUNNING** 7436592 (mrs 0.3/0.5 + PMC + domains) |
| **7434157** | `covid_v7_pmc_hotspot` | **TIMEOUT** (24h; best kept) | ep36 / 1.7316 | **RUNNING** 7436771 (mrs 0.3/0.5 + PMC + domains) |
| 7431995/96 | Track A v5/v6 | **RUNNING** (~36.5h / 48h) | — | No JSON; still K=100 group1; Coverage@K **—** |

Artifacts: `results/wave5_pull/` (logs + ckpt meta + existing v5_mutrec JSONs). Primary Table 5 draft row remains **covid_v5_mutrec @ mrs=0.5** until deepmut/h3n2/v7 JSONs land.

---

## TEMP paste-ready numbers (2026-08-06 Wave-4 siteaa / J.1)

Sources: `table_empirical_main.csv` (job 7350020 → `make_table.py`), [`WAVE4_EVAL_RESULTS.md`](../checkpoints/WAVE4_EVAL_RESULTS.md), `results/wave4_pull/` (+ Wave-3 temporal/clade still in `checkpoints/wave3_pull/`). `covid_v5_deepmut` + `covid_v7_pmc_hotspot` trained (TIMEOUT, best kept); enrichments **pending** 7436592 / 7436771 / 7436591 — do not paste until JSON exists.

### Table 2 — H3N2 tree generation (**TEMP / TOPO-OK**)

From `benchmarks/results/tables/table_empirical_main.csv` (n_roots=97). Tree-KL / Split-KL are `nan` → **—**. Adapted rows = pool-adapted only ([`BLOCKERS.md`](BLOCKERS.md)).

| Method (CSV / paper) | RF ↓ | Quartet ↓ | Branch W1 ↓ | Terminal edit ↓ | Notes |
|----------------------|-----:|----------:|------------:|----------------:|-------|
| NeutralBD / Neutral CTMC+BD | **0.9824 ± 0.0015** | **0.7110 ± 0.0067** | **7.22e−4 ± 6.45e−5** | **0.0148 ± 0.0048** | best RF |
| EmpiricalBD / JTT+BD proxy | **0.9839 ± 0.0014** | **0.7119 ± 0.0066** | **7.22e−4 ± 6.45e−5** | **0.0151 ± 0.0048** | not literal JTT CTMC — caption carefully |
| PLMPrior / pLM prior | **0.9823 ± 0.0014** | **0.7128 ± 0.0066** | **7.22e−4 ± 6.45e−5** | **0.0151 ± 0.0048** | |
| artreeformer_adapted | **0.9891 ± 0.0013** | **0.7126 ± 0.0065** | **8.61e−4 ± 7.10e−5** | **0.0151 ± 0.0048** | adapted pool |
| phylovae_adapted | **0.9900 ± 0.0012** | **0.7160 ± 0.0061** | **8.12e−4 ± 6.73e−5** | **0.0151 ± 0.0048** | adapted pool |
| **TreeSBM** (`treesbm`) | **0.9841 ± 0.0017** | **0.7179 ± 0.0061** | **8.64e−4 ± 7.65e−5** | **0.0136 ± 0.0048** | best Quartet among filled; TE slightly best |
| PhylaFlow / DCA (`phylaflow` native) | — | — | — | — | wired; needs Betty clone + **H3N2** bank (not DS1–8) + pools (`_paste_betty_phylaflow_native.sh`) |

### Table 5 — Viral mut / antigenic (**TEMP**; primary draft row = v5 mrs=0.5)

COVID geo enrichment, n_trees=20. Siteaa re-eval jobs **7431993/7431994** (`results/wave4_pull/`). Primary **TEMP** TreeSBM cell: **mut_recovery ≈ 0.0690** @ mrs=0.5 (Wave-3 was 0.0796 — re-seed drift). Cons ≪ 0.98; EVEscape still fails vs GT.

| Ckpt | mrs | mut (**TEMP**) | cons | id | site_recall | aa_acc\|hit | site×aa | EVEscape m/GT/rand |
|------|----:|---------------:|-----:|---:|------------:|------------:|--------:|--------------------|
| covid_v3_cons | 0.3 | 0.0201 | 0.898 | 0.688 | — | — | — | −0.005 / 0.047 / ~0 |
| covid_v3_cons | 0.5 | 0.0245 | 0.829 | 0.638 | — | — | — | +0.012 / 0.047 / ~0 |
| covid_v4_mutrec | 0.3 | 0.0428 | 0.899 | 0.691 | — | — | — | −0.019 / 0.047 / ~0 |
| covid_v4_mutrec | 0.5 | 0.0633 | 0.840 | 0.654 | — | — | — | −0.015 / 0.047 / ~0 |
| covid_v5_mutrec | 0.3 | 0.0402 | 0.904 | 0.701 | 0.186 | 0.273 | 0.051 | −0.011 / 0.047 / ~0 |
| **covid_v5_mutrec** | **0.5** | **0.0690** | 0.839 | 0.656 | **0.209** | **0.341** | **0.071** | −0.001 / 0.047 / ~0 |
| covid_v5_deepmut | 0.3 / 0.5 | — | — | — | — | — | — | **pending** job 7436592 |
| covid_v6_hotspot | 0.3 | 0.0554 | 0.906 | 0.704 | 0.183 | 0.336 | 0.061 | +0.015 / 0.047 / ~0 |
| covid_v6_hotspot | 0.5 | 0.0473 | 0.848 | 0.663 | 0.237 | 0.196 | 0.047 | +0.006 / 0.047 / ~0 |
| covid_v7_pmc_hotspot | 0.3 / 0.5 | — | — | — | — | — | — | **pending** job 7436771 (PMC + domains) |
| h3n2_v3_lit_hotspot | 0.3 / 0.5 | — | — | — | — | — | — | **pending** 7436591 (flu lit frac) |

Paper columns Coverage@1000 / Future rank / Min dist / true antigenic recall: **—** (Track A not finished). Influenza TreeSBM mut proxy: H1N1 temporal mut=**0.0856** (below); H3N2 lit-hotspot enrichment pending. HIV: **—**.

### Table 3 — Future lineage (**TEMP** mut proxies; H3N2 curves filled)

**H3N2 Coverage@K curves (job 7439926 COMPLETED):** `benchmarks/results/coverage_curves_h3n2_N16.csv` (N=16, 5 roots, ε=0.02, K=10..100). Plots: `benchmarks/results/plots/coverage_curves_h3n2_N16.{png,pdf}`. At K=100 all methods hit coverage=1.0; site_recall leader **artreeformer_adapted 0.848**, then TreeSBM **0.740** / NeutralBD **0.732** / PLMPrior **0.553**; best mean_min_edit **TreeSBM 2.35**.

COVID Track A Coverage@K still separate — jobs **7431995/7431996** **TIMEOUT**ed without full JSON set (see `results/wave5_pull/track_a/`).

| Split / ckpt | Protocol | n | mut_recovery | cons | identity | Coverage@K |
|--------------|----------|--:|-------------:|-----:|---------:|------------|
| H1N1 leaf-holdout | leaf-holdout | 284 | 0.0426 | 0.718 | **0.0808†** | **—** |
| └ test | leaf-holdout | 28 | 0.0407 | 0.722 | 0.0805† | **—** |
| **h1n1_temporal_v1** | enrichment mrs=0.3 | 20 | **0.0856** | 0.881 | **0.870** | **—** |
| **covid_temporal_v1** | enrichment mrs=0.3 | 16 | 0.0265 | 0.899 | 0.626 | **—** (EVEscape model≈GT +0.040/0.039) |
| covid_cladeholdout | leaf-holdout | 95 | 0.0485 | 0.772 | **0.097†** | **—** |
| covid_v5 Track A | `track_a.py` K∈{100,500,1000} | 20 | — | — | — | **— (pending)** |
| covid_v6 Track A | same | 20 | — | — | — | **— (pending)** |

† Suspect identity vs enrichment ~0.6–0.87 — draft on mut/cons only until protocol debug.

### Appendix J.1 — GraphTF probes (**FINAL**)

Artifacts: `results/wave4_pull/transformer_val_h3n2_j1/` (job **7434252**, `h3n2_v2`, n=222) and `results/wave4_pull/transformer_val_covid_v5_j1/` (job **7431998**, `covid_v5_mutrec`, n=335). Primary columns: plm / esm_br / graph_tf (full wide tables in artifacts).

**H3N2**

| Task | metric | plm | esm_br | graph_tf | topo | node_enc | graph_rand |
|------|--------|----:|-------:|---------:|-----:|---------:|-----------:|
| mutated_position | acc | **0.695** | 0.665 | 0.062 | 0.037 | 0.184 | 0.186 |
| ancestral_aa_identity | acc | **0.800** | 0.798 | 0.188 | 0.131 | 0.399 | 0.388 |
| leaf_vs_internal | auroc | 0.540 | **1.00** | 0.996 | 1.00 | 1.00 | 1.00 |
| subtree_size | r2† | −0.019 | 0.107 | **0.499** | 0.090 | 0.140 | 0.319 |
| parent_child_pair | auroc | 0.983 | 1.000 | **1.000** | 1.000 | 0.999 | 0.999 |
| num_children | acc | 0.573 | 0.941 | 0.927 | **0.941** | 0.941 | 0.941 |

**COVID**

| Task | metric | plm | esm_br | graph_tf | topo | node_enc | graph_rand |
|------|--------|----:|-------:|---------:|-----:|---------:|-----------:|
| mutated_position | acc | **0.506** | 0.444 | 0.010 | 0.032 | 0.083 | 0.013 |
| ancestral_aa_identity | acc | **0.777** | 0.772 | 0.129 | 0.130 | 0.341 | 0.310 |
| leaf_vs_internal | auroc | 0.543 | **1.00** | 1.00 | 1.00 | 1.00 | 1.00 |
| subtree_size | r2† | −0.043 | 0.084 | **0.398** | 0.078 | 0.178 | 0.363 |
| parent_child_pair | auroc | 0.929 | **1.000** | 0.999 | 0.999 | 0.996 | 0.999 |
| num_children | acc | 0.532 | 0.972 | 0.974 | **0.976** | 0.976 | 0.975 |

† Primary metric code prefers spearman; these runs logged r2 only. Sequence probes: plm ≈ esm_br ≫ graph_tf. Topology: graph_tf leads subtree_size; leaf/parent nearly saturated for all non-plm.

---

## Table → script → metric map

| Paper table | What it reports | Source command | Artifact(s) | Key metrics |
|-------------|-----------------|----------------|-------------|-------------|
| **Table 2** (H3N2 tree gen) | Topology + sequence quality of generated trees vs held-out roots | `sbatch --qos=mig-max scripts/slurm_baselines.sh [ckpt]` → `benchmarks/make_table.py` | `benchmarks/results/results_baselines.csv` → `benchmarks/results/tables/table_empirical_main.csv` (+ `.tex`); sim track if run | RF, Quartet (`tqdist`), branch-length / sequence distances per `run_table.py`. Rows: NeutralBD, EmpiricalBD, PLMPrior, `*_adapted`, TreeSBM |
| **Table 3** (future lineage) | Recovery of never-seen / future leaves | Leaf-holdout: `sbatch scripts/slurm_h1n1_leafholdout_eval.sh`. Temporal: train/eval on `data/h1n1_temporal`, `data/h3n2` or `data/h3n2_temporal_forecast`, `data/covid_temporal` ([`SPLITS.md`](SPLITS.md)). COVID clade: `data/covid_cladeholdout` + `eval_leaf_holdout.py --max-seq-len 1280`. **H3N2 Coverage@K curves (locked):** `sbatch scripts/slurm_coverage_curves.sh` — 5 held-out roots, ε=0.02, `site_recall`, K=10..100, methods NeutralBD / PLMPrior / artreeformer_adapted / TreeSBM (`h3n2_v3_lit_hotspot`); plan [`TABLE3_BASELINE_COVERAGE_PLAN.md`](../checkpoints/TABLE3_BASELINE_COVERAGE_PLAN.md). COVID Track A: `python benchmarks/track_a.py …` | `benchmarks/results/coverage_curves_h3n2_N16.csv` (+ R plot); Track A JSON under `benchmarks/results/`; leaf-holdout JSON | `site_recall`, Coverage@ε, mean_min_edit vs K; mut_recovery / cons / identity proxies |
| **Table 4** | — | **Skip entirely** | — | — |
| **Table 5** (viral mut / antigenic) | COVID spike mut recovery + escape enrichment (+ EVE) | `sbatch --qos=mig-max scripts/slurm_eval_covid.sh` (v3); v4: `scripts/slurm_eval_covid_v4.sh`. mrs grid: `scripts/slurm_inference_sweep.sh` / `slurm_inf_sweep.sh`. EVE: `sbatch scripts/slurm_eval_eve.sh <ckpt> <data> <eve.pt> [L]` | `checkpoints/eval_enrichment_covid_*.json`; `checkpoints/inference_sweep_*/eval_mrs*.json`; `checkpoints/eval_eve_*.json` | `mut_recovery`, `cons_retention`, best_match_id, EVEscape (model / GT / rand); EVE recovered-mut vs random |
| **Table 6** (Ab affinity) | CDR mut / SHM / lineage RF | `scripts/prepare_ab_fasta.py` → IgBLAST → trees (see ICLR §B) | `data/ab_prep/` (shards); trees **not built yet** | all **—** until clone trees exist |

Table 7–8 (later): arch loop in [`MUT_RECOVERY_PLAYBOOK.md`](MUT_RECOVERY_PLAYBOOK.md) + S+H26 ideas in [`MUT_RECOVERY_SH26_ITERATIONS.md`](MUT_RECOVERY_SH26_ITERATIONS.md), then re-run Table 2/5.

---

## Current number status (Wave-5 pull, 2026-08-07)

From [`WAVE5_EVAL_RESULTS.md`](../checkpoints/WAVE5_EVAL_RESULTS.md) (`results/wave5_pull/`) + Wave-4 siteaa/J.1; temporal/clade still Wave-3. Paste-ready copy is in **TEMP paste-ready numbers** above.

| Table | Cell / row | Status | Real numbers (artifact) |
|-------|------------|--------|-------------------------|
| **T2** H3N2 empirical means (CSV) | BD / PLM / TreeSBM / adapted | **TEMP pasted** | n_roots=97; see TEMP Table 2 |
| **T2** PhylaFlow native (`phylaflow`) | — | **Empty** | Harness wired; train on H3N2 train (not DS1–8) — Betty: `bash scripts/_paste_betty_phylaflow_native.sh` |
| **T3** temporal / clade mut proxies | mut / cons / identity | **TEMP pasted** | Coverage@K still **—** |
| **T3** Coverage@K | Track A v5/v6 | **Empty / pending** | 7431995/96 RUNNING ~36.5h; no `track_a_*_K*.json`; ETA none before wall |
| **T3** H3N2 coverage curves | `coverage_curves.py` (5 roots, ε=0.02, site_recall, +TreeSBM) | **DONE** job **7439926** | CSV `coverage_curves_h3n2_N16.csv` + plots; @K=100: cov=1.0 all; site_r ART=**0.848** / TreeSBM=**0.740** / Neutral=**0.732** / PLM=**0.553**; edit TreeSBM best **2.35** |
| **T5** covid_v5 mrs=0.5 + site×AA | mut≈0.069 | **TEMP pasted** | mut=**0.0690** site_r=0.209 aa\|hit=0.341 |
| **T5** covid_v5_deepmut / h3n2 lit | enrichment | **pending evals** | jobs **7436592** / **7436591** |
| **T5** covid_v7 PMC hotspot | train TIMEOUT; eval | **RUNNING** 7436771 | best ep36 val=1.7316; mrs 0.3/0.5 |
| **T5** EVE Marks | — | **Empty** | Needs `prepare_eve_scores.py` → `.pt` then `slurm_eval_eve.sh` |
| **J.1** GraphTF H3N2+COVID | probes | **FINAL** | `transformer_val_h3n2_j1/` (7434252) + `transformer_val_covid_v5_j1/` (7431998) |
| **T6** Antibody | all | **Empty** | Ab prep/IgBLAST scaffolded; no trees yet |

Candidate Table 5 rows: prefer **v5 @ mrs=0.5** for max geo enrichment mut (**0.069** after siteaa re-eval; Wave-3 was 0.080) — cons≪0.98. Compare to deepmut/h3n2 lit/v7 once 7436591/92/7436771 finish.

---

## How to fill after longer baseline / eval runs

### Table 2 (baselines, 24h)

```bash
# On Betty (login shell so SLURM is on PATH):
cd ~/DiscreteTreeFlows
# Sync latest scripts from laptop if needed, then:

# H3N2 + TreeSBM (default paths in script):
sbatch --qos=mig-max scripts/slurm_baselines.sh checkpoints/h3n2_v2/best.pt

# Or BD/PLM only:
sbatch --qos=mig-max scripts/slurm_baselines.sh
```

Produces `benchmarks/results/results_baselines.csv` and aggregates under `benchmarks/results/tables/`.  
If Quartet is NaN: install/verify `tqdist` (script preamble tries); leave “—” in the paper, never fabricate.  
Adapted rows appear only when pools exist under `benchmarks/external_pools/sampled/` ([`EXTERNAL.md`](EXTERNAL.md), [`BLOCKERS.md`](BLOCKERS.md)).

**Resubmit after timeout:** same `sbatch` line; raise `--time` in the script header (currently `24:00:00`) or override:

```bash
sbatch --qos=mig-max --time=24:00:00 scripts/slurm_baselines.sh checkpoints/h3n2_v2/best.pt
```

### Table 3

```bash
sbatch scripts/slurm_h1n1_leafholdout_eval.sh
# After temporal/clade data exists:
python scripts/eval_leaf_holdout.py --data data/covid_cladeholdout \
  --checkpoint checkpoints/covid_v4_mutrec/best.pt --max-seq-len 1280
python benchmarks/track_a.py --checkpoint checkpoints/covid_v4_mutrec/best.pt \
  --data data/covid/test --K 10 --n-groups 20 \
  --out benchmarks/results/track_a_covid_v4.json
```

### Table 5

```bash
sbatch --qos=mig-max scripts/slurm_eval_covid.sh
sbatch --qos=mig-max scripts/slurm_eval_covid_v4.sh   # when v4 ckpt exists
sbatch --qos=mig-max scripts/slurm_inference_sweep.sh
# EVE (after prepare_eve_scores.py):
sbatch --qos=mig-max scripts/slurm_eval_eve.sh <ckpt> <data> <eve.pt> 1280
```

Pull JSONs locally; only then paste means into the paper. Confirm `mutation_rate_scale`, `n_steps`, `max_trees`, checkpoint name, and job ID in each artifact `summary`.

### Table 6 (antibody long pole — prep only)

```bash
# On Betty after syncing scripts + data/Homo_sapiens.fasta:
cd ~/DiscreteTreeFlows
mkdir -p logs data/ab_prep
# Parse + shard + filter (CPU; ~1–2h for full 3.4M):
sbatch scripts/slurm_ab_prepare_fasta.sh
# Or dry-run locally first:
python scripts/prepare_ab_fasta.py --input data/Homo_sapiens.fasta \
  --out-dir data/ab_prep --max-seqs 5000 --n-shards 4

# After shards exist + IgBLAST/IMGT refs installed on Betty:
sbatch scripts/slurm_ab_igblast_clones.sh   # scaffold; edit IGBLAST_BIN/GERMLINE_DIR
# Trees / train / T6 metrics: NOT done — see ICLR_TABLE_FILL_PLAN §B
```

---

## Aggregation helpers

- Long → paper CSV/LaTeX: `benchmarks/make_table.py --results … --out … --ci se|boot`
- Track A JSON: `benchmarks/track_a.py`
- EVEscape / mut-cons JSON: `scripts/eval_evescape_enrichment.py`
- Leaf-holdout JSON: `scripts/eval_leaf_holdout.py`
- EVE (Marks) JSON: `scripts/eval_eve_baseline.py`
- Ab FASTA prep: `scripts/prepare_ab_fasta.py`
- Ab CDR/SHM metric stubs: `scripts/ab_t6_metrics.py` (stubs only)

## Checklist before pasting into the paper

1. Confirm checkpoint name + commit / job ID in the artifact `summary`.
2. Confirm eval flags (`mutation_rate_scale`, `n_steps`, `max_trees`) match the caption.
3. Confirm split protocol (geo vs temporal vs clade) matches [`SPLITS.md`](SPLITS.md) and the caption.
4. Adapted rows only if pools exist under `benchmarks/external_pools/sampled/` — never claim native forward gen (`BLOCKERS.md`).
5. Quartet NaN → “—” + footnote missing `tqdist`; never fabricate.
6. Mark mut_recovery cells **TEMP** until DEEP_MUT / better recipe (not queued in this pass).
