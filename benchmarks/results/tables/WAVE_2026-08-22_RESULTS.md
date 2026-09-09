# Retrain + VaxSeer + Ab SHM wave — results (2026-08-22)

Pulled from Betty sacct + artifacts. **No invented numbers.** Paper ckpts were not overwritten.

Snapshot **2026-08-23**: COVID mutlin trains **COMPLETED**; mutlin enrich **7794867–7794874** COMPLETED; Ab Recipe B **7799979–7799981** COMPLETED. Paper ckpts were not overwritten.

---

## Job board (IDs from 2026-08-20 submit)

| ID | Name | State | Elapsed | Notes |
|---:|---|---|---|---|
| 7720891–902 | enrich mrs 0.3/0.5/1.0 × 4 viruses | **COMPLETED** (JSONs on disk) | ~20–90 min | `checkpoints/viral_eval_sweep/` |
| 7720903 | covid Cov K=100 | COMPLETED | 6:05 | |
| 7720904 | covid Cov K=500 | **TIMEOUT** 12h | empty CSV | |
| 7720905 | h1n1 Cov K=100 | COMPLETED | 2:49 | |
| 7720906 | h1n1 Cov K=500 | **TIMEOUT** 12h | empty CSV | |
| 7720907 | hiv_geo Cov K=100 | COMPLETED | 1:33 | |
| 7720908 | hiv_geo Cov K=500 | COMPLETED | 7:16 | only K=500 that finished |
| 7720909 | h3n2 Cov K=100 | COMPLETED | 4:21 | |
| 7720910 | h3n2 Cov K=500 | **TIMEOUT** 12h | empty CSV | |
| 7720911 | covid_v6_mutlin_b0 | **COMPLETED** | 2-08:53 | `best.pt` epoch 67 val **0.651** |
| 7720912 | h1n1_v3_mutlin_b0 | COMPLETED | 20:02 | ep 74 val **3.551** |
| 7720913 | h3n2_v4_mutlin_b0 | COMPLETED | 10:36 | ep 89 val **5.151** |
| 7720914 | hiv_geo_v2_mutlin_b0 | COMPLETED | 4:12 | ep 14 val **76.67** (high vs flu; inspect before using) |
| 7720915 | covid_v6_mutlin_b025 | **COMPLETED** | 2-07:09 | `best.pt` epoch 67 val **0.810** |
| 7720916 | h1n1_v3_mutlin_b025 | COMPLETED | 20:24 | ep 74 val **3.777** |
| 7720917 | h3n2_v4_mutlin_b025 | COMPLETED | 10:39 | ep 57 val **6.087** |
| 7720918 | hiv_geo_v2_mutlin_b025 | COMPLETED | 4:03 | ep 14 val **79.12** |
| 7720919 | ab_oas_v2_shm train | COMPLETED | 6:39 | ep 56 val **87.78**; β=0, shm boost=2, fwr=0.5, AID on |
| 7720920 | ab_oas_v2 Rod.82 | COMPLETED | 3:18 | 1640 samples `treesbm_ab_oas_v2` |
| 7720921 | h3n2_treeviz | COMPLETED | 0:47 | groups **55** and **40** matched |
| 7720923 | ab entropy audit | COMPLETED | 0:00:16 | H_CDR/H_FWR OAS **1.627**, Rod **2.024** |
| 7720924 | h3_varpath | COMPLETED | 0:00:01 | wrote PENDING (filename mismatch; see below) |

---

## 1. Eval-only mrs sweep (current paper ckpts)

Free-gen enrichment (`eval_evescape_enrichment.py`, n_trees=20 except HIV 14). Primary: `mut_recovery`, `aa_acc_given_hit` (= aa\|hit), `cons_retention`. **These JSONs have no `clade_recall` / `mean_min_edit`** — those come from coverage CSVs (Table 2).

| virus | ckpt | mrs | mut_recovery | aa\|hit | site_recall | cons_retention | pLM NLL |
|---|---|---:|---:|---:|---:|---:|---:|
| COVID | covid_v5_mutrec | 0.3 | 0.0446 | 0.290 | 0.173 | 0.904 | 0.457 |
| COVID | covid_v5_mutrec | **0.5** | **0.0690** | **0.341** | 0.209 | 0.839 | 0.457 |
| COVID | covid_v5_mutrec | 1.0 | 0.121 | 0.233 | 0.425 | 0.711 | 0.452 |
| H1N1 | h1n1_v2_lit_hotspot | 0.3 | 0.109 | 0.626 | 0.164 | 0.971 | 0.462 |
| H1N1 | h1n1_v2_lit_hotspot | **0.5** | **0.154** | **0.618** | 0.241 | 0.952 | 0.461 |
| H1N1 | h1n1_v2_lit_hotspot | 1.0 | 0.174 | 0.430 | 0.355 | 0.907 | 0.461 |
| H3N2 | h3n2_v3_lit_hotspot | 0.3 | 0.124 | **0.669** | 0.199 | 0.953 | 0.461 |
| H3N2 | h3n2_v3_lit_hotspot | 0.5 | 0.100 | 0.439 | 0.253 | 0.918 | 0.463 |
| H3N2 | h3n2_v3_lit_hotspot | 1.0 | **0.232** | 0.559 | 0.413 | 0.849 | 0.465 |
| HIV geo | hiv_geo_v1 | 0.3 | 0.138 | 0.404 | 0.258 | 0.946 | 0.501 |
| HIV geo | hiv_geo_v1 | **0.5** | **0.154** | 0.310 | 0.386 | 0.896 | 0.503 |
| HIV geo | hiv_geo_v1 | 1.0 | 0.197 | 0.324 | 0.525 | 0.785 | 0.499 |

**Read:** raising mrs trades conservation for mut/site recall. COVID aa\|hit peaks at mrs=0.5 then drops at 1.0. H3 mut peaks at mrs=1.0 but aa\|hit is best at 0.3. Do not bake mrs=1.0 into train without checking aa\|hit.

---

## 2. Coverage K / ε (same paper ckpts, 5 locked roots, N=16)

`coverage_obs_e*` = fraction of **observed** leaves within Hamming ε of some generated leaf. K=500 timed out for COVID/H1/H3 (12h wall); HIV geo K=500 finished.

| virus | K | Cov@ε=1 | Cov@ε=2 | Cov@ε=3 | unique@ε2 | mean_min_edit | mut_recovery | aa\|hit | clade_recall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| COVID | 100 | 0.775 | **0.775** | 0.788 | 0.750 | 1.41 | 0.0008 | 1.00 | **0.000** |
| H1N1 | 100 | 0.688 | **0.800** | 0.925 | 0.761 | 1.25 | 0.0565 | 1.00 | **0.175** |
| H3N2 | 100 | 0.625 | **0.775** | 0.800 | 0.749 | 2.11 | 0.266 | 1.00 | **0.300** |
| HIV geo | 100 | 0.138 | **0.188** | 0.225 | 0.188 | 63.09 | 0.0262 | 0.614 | **0.078** |
| HIV geo | 500 | 0.150 | **0.188** | 0.238 | 0.188 | 62.89 | 0.0579 | 0.742 | **0.144** |

COVID Cov@ε=2@K=100 stays **0.775** (same plateau as `COVERAGE_WHAT_MOVES_IT.md`). HIV Cov@ε=2 does **not** move from K=100→500; clade_recall does (0.078→0.144).

---

## 3. Viral retrain (`*_mutlin`, β ∈ {0, 0.25} site_local)

New dirs only. Enrichment mrs=0.5: jobs **7794867–7794874** COMPLETED. JSONs in `checkpoints/viral_eval_sweep/*mutlin*_enrich.json`; compact means: [`viral_mutlin_enrich_summary.json`](viral_mutlin_enrich_summary.json). **No `clade_recall` in these files** (same as §1 enrich). Compare to paper-ckpt **mrs=0.5** rows in §1.

| ckpt | β | epoch | val_loss | n | mut_recovery | aa\|hit | site_recall | cons_retention | vs paper mrs=0.5 mut / aa\|hit |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| covid_v6_mutlin_b0 | 0 | 67 | 0.651 | 20 | **0.0698** | 0.277 | 0.219 | 0.835 | paper 0.0690 / **0.341** — mut tied, aa\|hit worse |
| covid_v6_mutlin_b025 | 0.25 | 67 | 0.810 | 20 | 0.0554 | 0.306 | 0.206 | 0.853 | mut down vs paper |
| h1n1_v3_mutlin_b0 | 0 | 74 | 3.551 | 20 | 0.124 | 0.558 | 0.184 | 0.953 | paper **0.154 / 0.618** — both worse |
| h1n1_v3_mutlin_b025 | 0.25 | 74 | 3.777 | 20 | 0.139 | **0.628** | 0.196 | 0.969 | mut still below paper; aa\|hit ≈ paper |
| h3n2_v4_mutlin_b0 | 0 | 89 | 5.151 | 20 | 0.118 | 0.416 | 0.250 | 0.919 | paper 0.100 / 0.439 — mut up, aa\|hit similar |
| h3n2_v4_mutlin_b025 | 0.25 | 57 | 6.087 | 20 | **0.182** | **0.685** | 0.252 | 0.931 | best mutlin pair here; beats paper mrs=0.5 (not mrs=1.0 mut 0.232) |
| hiv_geo_v2_mutlin_b0 | 0 | 14 | 76.67 | 14 | 0.166 | 0.359 | 0.357 | 0.890 | paper 0.154 / 0.310 — slight mut/aa up; val still huge |
| hiv_geo_v2_mutlin_b025 | 0.25 | 14 | 79.12 | 14 | 0.149 | 0.348 | 0.326 | 0.922 | ≈ paper mut |

β=0.25 still **worsens val** on every pair. On free-gen enrich it **helps H1/H3 aa\|hit and H3 mut**, **hurts COVID mut**. HIV val ≫ flu — do not treat hiv_geo_v2 as a paper replacement.

---

## 4. Flu named-clade figure

Tree viz **7720921** completed: closest matched groups **55** then **40** (combo 0.545 / 0.547; RF=1; term-edit ~0.125; mean best-id ~0.90). Files: `results/h3n2_tree_viz/group_{055,040}_generated_matched.fasta`.

Follow-up **7720924** finished in **1s** with `h3n2_variant_path_PENDING.md`: the path script looked for `*_gen_leaves.fasta` / `*_matched_leaves.fasta`, not `*_generated_matched.fasta`. Re-run `build_flu_variant_path_figure.py` against the viz filenames (no new GPU job needed).

---

## 5. Antibody SHM

### Entropy audit (7720923)

| set | n_seqs | n_CDR / n_FWR | H_CDR | H_FWR | H_CDR/H_FWR |
|---|---:|---|---:|---:|---:|
| OAS train tips (80 groups) | 1873 | 27 / 133 | 2.245 | 1.380 | **1.627** |
| Rod.82 leaves | 964 | 27 / 133 | 2.054 | 1.015 | **2.024** |

CDR is more diverse than FWR under the IMGT mask — supports Recipe A stay boost on CDR∪AID, not ESM β>0.

### OAS v2 train + Rod.82

- Ckpt: `checkpoints/ab_oas_1m_v2_shm/best.pt` (does not overwrite `ab_oas_1m_v1`)
- Samples: `antibody_benchmark/results/samples/treesbm_ab_oas_v2` (1640 files, 82 families × 20)

| model | CDR recall | SHM load err | Cov@ε2 | Cov@ε5 | term-div err | Root→leaf W1 | site-freq ρ |
|---|---:|---:|---:|---:|---:|---:|---:|
| treesbm_ab_oas **v1** (N=20, paper) | 0.150 | 0.077 | 0.024 | — | 9.89 | — | — |
| treesbm_ab_oas v1 (re-eval 7720920, larger pool†) | 0.450 | 0.078 | 0.021 | 0.087 | 9.88 | 9.50 | 0.020 |
| **treesbm_ab_oas_v2** (N=20, this wave) | **0.162** | **0.078** | **0.021** | **0.087** | **9.98** | 9.53 | 0.040 |
| Neutral / Thrifty / CoSiNE (same T6 pass) | 0.75 / 0.75 / 0.66 | 0.037 / 0.028 / 0.037 | ~0.018–0.021 | 0.058 / 0.085 / 0.078 | 13.0 / 2.71 / 3.02 | — | — |

† 7720920’s Table 6 pass scored v1 with `n_gen_leaves=400` / 100 rollouts on some families, so CDR 0.450 is **not** comparable to the paper N=20 row. Fair v1 vs v2 at N=20: **0.150 → 0.162**. Recipe A did not reach the 0.4 target.

### OAS v3 Thrifty Q0 (Recipe B) — 7799979 / 7799980 / 7799981 COMPLETED

- Ckpt: `checkpoints/ab_oas_1m_v3_thrifty/best.pt` (epoch 56, val **71.33**, `r0_backend=thrifty_aa`, β=0, no Recipe A stay prior). Does not overwrite v1/v2 or ESM `group_*_ref_rates.pt`.
- Caches: 397 `group_*_ref_rates_thrifty.pt` under `data/ab_clones_1m/train`.
- Samples: `antibody_benchmark/results/samples/treesbm_ab_oas_v3` (1640 files, 82×20).
- Fair T6 N=20: [`table6_ab_oas_v2v3.md`](table6_ab_oas_v2v3.md). Primary rollout: Root→leaf W1 **10.13**, site-freq ρ **0.035**.

| model | CDR recall | SHM load err | Cov@ε2 | Cov@ε5 | term-div err |
|---|---:|---:|---:|---:|---:|
| treesbm_ab_oas **v1** (N=20, paper) | 0.150 | 0.077 | 0.024 | — | 9.89 |
| treesbm_ab_oas_v2 Recipe A (N=20) | 0.162 | 0.078 | 0.021 | 0.087 | 9.98 |
| **treesbm_ab_oas_v3** Recipe B (N=20) | **0.175** | **0.083** | **0.018** | **0.087** | **11.13** |
| Neutral / Thrifty / CoSiNE | 0.69 / 0.75 / 0.66 | 0.054 / 0.037 / 0.037 | ~0.015–0.021 | 0.059 / 0.058 / 0.076 | 17.4 / 13.0 / 3.02 |

Recipe B is a small CDR bump over A (**0.162 → 0.175**) and still far from Thrifty **0.75**. SHM-load and terminal-diversity errors got slightly worse. Did **not** hit the 0.4 CDR target.

### OAS v4 cosine RateHeads (Recipe C) — train **7816338** COMPLETED, Rod.82 **7816339** COMPLETED

Site-local CNN mut-head, parent→child CE, R0 mix init σ(−4), ESM-2 caches reused, β=0. Ckpt `checkpoints/ab_oas_1m_v4_cosinehead/best.pt` (epoch 47, val **134.73** — much worse than v3 val 71.33). Samples `treesbm_ab_oas_v4` (1640 files). Fair T6: [`table6_ab_oas_v2v3v4.md`](table6_ab_oas_v2v3v4.md). Primary: Root→leaf W1 **6.62**, site-freq ρ **−0.057**. The first Rod.82 T6 pass omitted the v4 row; this table is a re-score of the same samples (does not overwrite v1–v3 sample dirs).

| model | CDR recall | SHM load err | Cov@ε2 | Cov@ε5 | term-div err |
|---|---:|---:|---:|---:|---:|
| treesbm_ab_oas **v1** (N=20, paper) | 0.150 | 0.077 | 0.024 | — | 9.89 |
| treesbm_ab_oas_v2 Recipe A (N=20) | 0.162 | 0.078 | 0.021 | 0.087 | 9.98 |
| treesbm_ab_oas_v3 Recipe B (N=20) | **0.175** | 0.083 | 0.018 | 0.087 | 11.13 |
| **treesbm_ab_oas_v4** Recipe C (N=20) | **0.158** | **0.054** | **0.018** | 0.087 | **4.89** |
| Neutral / Thrifty / CoSiNE | 0.69 / 0.75 / 0.66 | 0.054 / 0.037 / 0.037 | ~0.015–0.021 | 0.059 / 0.058 / 0.076 | 17.4 / 13.0 / 3.02 |

Recipe C **does not** improve CDR vs B (0.175 → **0.158**, back near A). SHM-load error and terminal-diversity error improved (0.083→0.054, 11.13→4.89). Site-frequency ρ went **negative**. Still far from Thrifty CDR **0.75**.

---

## 6. VaxSeer blocker (why generated vs observed \(p_Y\) is empty)

This is **not** a Slurm failure. The jobs in `squeue` (eval + retrain) never included a VaxSeer GPU/CPU scorer, because **live dominance LM weights never downloaded**.

What we have:

| Piece | Status |
|---|---|
| Clone `$LABHOME/third_party/vaxseer` | done |
| Official **results tar** (`people.csail.mit.edu/wxsh/vaxseer/results.tar.gz`) | done — `dominance_prediction/lm/.../test_results.csv` |
| Adapter + `eval_vaxseer_dominance.py` | done |
| **lm checkpoints** (`runs/flu_lm/.../*.ckpt`) | **missing** |

Why weights are missing:

1. Upstream `download_models_from_dropbox.py` calls `input("Please input your dropbox access token: ...")`. Non-interactive SSH/sbatch **cannot** type that token. There is no env-var bypass in their script.
2. The Dropbox folder is [shared via a personal link](https://www.dropbox.com/scl/fo/7d94eqsii2h1jdm5l7mm6/h?rlkey=1n1wafyuapwx5a4c04jc0y7cs&dl=0); the downloader uses the Dropbox Python API (`sharing_get_shared_link_file`) **after** a user token.
3. Year grid is **2012–2021 only** (not 2024). Even after download, TreeSBM H3 temporal 2024 test is not year-aligned until you pick a VaxSeer year and regenerate from pre-Feb roots.
4. The extracted CSVs map **GISAID EPI ids → score**, not AA strings. They can rank *observed* strains that already have EPI ids. **Generated TreeSBM leaves have no EPI id**, so they cannot be scored from the CSV. Scoring \(p_Y(g)\) needs the frozen GPT2-time LM ckpt.

**Unblock (you, once):** on a login node with a Dropbox app token:

```bash
cd /vast/projects/pranam/lab/nnori/third_party/vaxseer
# conda env with `dropbox` (treesbm env lacked it)
python download_models_from_dropbox.py --task lm --year 2018 --subtype a_h3n2 a_h1n1 --output_dir runs
```

Then point `VAXSEER_ROOT` / `VAXSEER_RUNS` at those ckpts and run `scripts/eval_vaxseer_dominance.py` on H3/H1 generated FASTAs vs season observed HA.

## 7. Why aa|hit is not doing what we want

`aa|hit` = P(generated AA = GT | root≠GT and generated≠root). Site recall can rise while this falls: `mut_recovery = site_recall × aa|hit`.

**What the numbers already showed (viral mutlin, mrs=0.5):** COVID aa|hit **0.28–0.31**; HIV **0.35**; H1 **0.56–0.63**; H3 β=0.25 **0.685** vs β=0 **0.416**. Paper COVID mrs 0.5→1.0: site recall up, aa|hit **0.341→0.233**. mrs is a volume knob, not an identity knob.

**Polymerase-stall biology vs this metric.** Homopolymer stutter, hairpins/pseudoknots, G-quadruplexes, and palindromic/direct repeats mainly cause **indels, frameshifts, backtracking, and template switches**. TreeSBM’s sampler (`mutate_sequence_independent`) is **per-site AA substitution only** — no codon, no NT, no indels. Those stall features cannot show up as “wrong AA at a hit site” except indirectly (if a real indel is scored as a pile of substitutions). Putting RT/RdRp stall motifs into the **AA** RateHead will not, by itself, make aa|hit look like CoSiNE/Thrifty.

**COVID matched trees (groups 3/9/36/40), proxy reverse-translate of the root AA** (most-common codon, not authentic SARS-CoV-2 NT):

| stratum | GT mut sites | site recall | aa\|hit | chem\|hit (BLOSUM-close) |
|---|---:|---:|---:|---:|
| all | 76831 | 0.302 | **0.212** | 0.301 |
| no stall motif | 70312 | 0.304 | 0.209 | 0.302 |
| any stall | 6519 | 0.283 | 0.238 | 0.294 |
| homopolymer | 5228 | 0.300 | 0.197 | 0.257 |
| tandem repeat | 1455 | **0.207** | 0.272 | 0.329 |
| palindrome | 366 | 0.303 | 0.550 | 0.550 |

Stall windows are **not** where aa|hit collapses (homopolymer 0.20 vs background 0.21). Tandem repeats are **under-hit** (site recall 0.21 vs 0.30). Palindrome aa|hit looks high but n is small. **chem|hit − aa|hit ≈ 0.09**: many hits are biochemically nearby AAs, not the observed residue. No G4 hits on this reverse-translated spike (expected: proxy NT ≠ genome NT).

**What is actually broken**

1. Destination mass is independent-site softmax over 20 AAs (residual ESM or Recipe C CNN). It does not condition on codon, neighboring NT, or polymerase processivity.
2. Scaling mrs / mutlin raises P(leave root) without sharpening P(AA | leave). Extra hits sample the same uninformed dest.
3. Recipe C trained parent→child CE at val **134.7** (v3 was **71.3**) and **lowered** Rod.82 CDR (0.175→0.158). A site CNN without NT context did not learn SHM or viral dest identity.
4. Authentic stall priors need **genome NT** (or a codon model that emits indels). Preferred-codon reverse translation is only a motif proxy.

**If we use that biology next:** NT-conditioned dest (codon + homopolymer/G4/repeat flags on real CDS), and/or an indel channel. Do not expect aa|hit to jump from mixing stall scores into the current AA Gillespie.

---

## What is still open

1. Mutlin enrich is in (mrs=0.5 only). Coverage K/ε and clade_recall were **not** rerun on mutlin ckpts.
2. Re-run flu path figure on `group_*_generated_matched.fasta`.
3. Dropbox token for VaxSeer lm; then generated-vs-observed dominance table.
4. Recipe C Rod.82 is in: CDR **0.158** (worse than B). aa|hit is a dest-identity problem, not a stall-motif site problem — see §7.
