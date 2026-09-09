# Antibody — three tracks (do not mix)

Pulled / status **2026-08-15** (Table 5 redo: Neutral SHM=JC69 ≠ Thrifty; AR tested; CoSiNE/TreeSBM filled; pLM Rod.82 in progress).

There are **three separate antibody experiments**. Metrics are not interchangeable.

---

## Track A — OAS 1M Homo_sapiens TreeSBM (… → **7627289** FAILED-with-ckpt / Rod.82 **7627919** COMPLETED)


**What was done**

1. Subsample **1,000,000** heavy-chain AA seqs from `data/Homo_sapiens.fasta` (3,435,378 read; 1,648,017 length/chain-filtered; 1,440 donors; 604,326 preferred-paired).
2. ANARCI IMGT: 999,362 / 1,000,000 OK CDR3 (`python_api`, 8 jobs).
3. Clones = `donor_id|v_family|cdr3`; min_size=16 max_size=64 → **611 kept clones** (633,164 raw keys; 295,675 no-donor dropped).
4. MAFFT + FastTree-LG + midpoint root; root seq = majority-leaf vote (**not germline ASR**).
5. **Donor-disjoint** clonal-lineage split (seed=42): train 397 / val 44 / test 170 groups.
6. CDR hotspot mask empirical from ANARCI: `results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt` — 27 hot sites / 160, frac=0.169, occupancy≥0.25.
7. **Postprocess fix (2026-08-13):** `scripts/ab_postprocess_clones_for_treesbm.py` wrote `*_rooted.nwk` + leaf+root `*_anc_aa.fasta` + `*_bl.json` (`numdate`=cum path length). Complete: train **397** / val **44** / test **170**.
8. Train: entropy ON + CDR mask + **`LAMBDA_BR=0`** → `checkpoints/ab_oas_1m_v1`.

| Step | Job | State | Elapsed | Notes |
|---|---:|---|---|---|
| preprocess | **7585306** | COMPLETED | 02:33:01 | genoa-std-mem; wrote `data/ab_t5_1m` + `data/ab_clones_1m` |
| train afterok | **7585307** | **FAILED** 1:0 | 00:02:59 | TreeDataset **0 trees** (missing rooted/bl) |
| postprocess + retrain | **7596655** | **FAILED** | 00:00:09 | pre-repair / bad launch |
| retrain after support-strip repair | **7618018** | **FAILED** | — | val=nan every epoch (gap internals); no `best.pt`; Rod.82 **7618019** DependencyNeverSatisfied |
| retrain after parent-fill fix | **7627289** | **FAILED** 1:0 | 05:15:26 | Wrote `best.pt` (val=**64.4053**, early-stop @ep122) then crashed on post-train PLM (test split had 0 `*_plm.pt`) |
| Rod.82 afterok | **7627290** | **CANCELLED** | — | `DependencyNeverSatisfied` — cancelled 2026-08-15 |
| precompute PLM/ref_rates | **7627915** | **COMPLETED** | 00:04:03 | all complete groups cached (train 397 / val 44 / test 170) |
| Rod.82 no-dep resubmit | **7627916** | **FAILED** | 00:00:53 | `aa_indices` missing for `use_mut_aa_emb` OAS ckpt |
| Rod.82 after treesbm fix | **7627919** | **COMPLETED** | 01:03:40 | `treesbm_ab_oas`, N=20, observed topo+BL → `table6_ab_track_c.*` |

See `results/ab_cdr_mask/OAS_1M_POSTPROCESS_RETRAIN.md`.

**Holdout (Track A only — not Rodriguez 82)**

| Split | n groups | n donors | Donors |
|---|---:|---:|---|
| train | 397 | 56 | PRJEB40825, PRJEB51634, PRJEB53053, PRJEB61178, SRR10596378, SRR10596382, SRR10596387, SRR10596395, SRR10596400, SRR10596404, SRR10596409, SRR12875362, SRR14978396, SRR14978397, SRR17729677, SRR17729679, SRR17729682, SRR17729689, SRR17729690, SRR17729691, SRR17729692, SRR17729694, SRR17729695, SRR17729699, SRR17729703, SRR17729706, SRR17729707, SRR17729708, SRR17729711, SRR17729714, SRR17729715, SRR17729716, SRR17729725, SRR17729726, SRR17778090, SRR17778131, SRR17778142, SRR17816506, SRR17816508, SRR17816511, SRR17816512, SRR17818204, SRR17818206, SRR17818209, SRR18356101, SRR18363731, SRR20210850, SRR20210852, SRR22865396, SRR24716320, SRR24716342, SRR24716354, SRR24716355, SRR24716367, SRR24800050, SRR25405110 |
| val | 44 | 7 | SRR17729673, SRR17778096, SRR17778097, SRR17778143, SRR17816507, SRR20210854, SRR24716366 |
| test (held out) | 170 | 10 | SRR10596391, SRR17729674, SRR17729678, SRR17729680, SRR17729700, SRR17729701, SRR17816505, SRR17816509, SRR17818205, SRR17818208 |

Donor overlap train∩val∩test = empty (donor-disjoint). Test is dominated by **SRR17818205** (124/170 clones); then SRR17816509 (14), SRR17816505 (12), SRR17818208 (10).

Test donor × n_clones: SRR17818205=124, SRR17816509=14, SRR17816505=12, SRR17818208=10, SRR17729674=2, SRR17729680=2, SRR17729700=2, SRR17729701=2, SRR10596391=1, SRR17729678=1.

Val donor × n_clones: SRR17778143=19, SRR17816507=18, SRR17729673=3, SRR17778096=1, SRR17778097=1, SRR20210854=1, SRR24716366=1.

**What is being tested (Track A):** TreeSBM on OAS *clonal phylogenies* (FastTree on clone MSAs), with a **CDR-weighted hotspot mask** and **site-entropy loss weighting**, on a donor-disjoint clone holdout. Ckpt on Betty: `checkpoints/ab_oas_1m_v1/best.pt`. Rod.82 eval **COMPLETED** as **7627919** (`treesbm_ab_oas`; Cov@e2=**0.0244**, SHM err=**0.0774**, CDR recall=**0.1496**, term-div err=**9.888**).

Paths: `data/Homo_sapiens.fasta` → `data/ab_t5_1m/` → `data/ab_clones_1m/{train,val,test}` · mask `results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt` · ckpt `checkpoints/ab_oas_1m_v1` · logs `logs/ab_oas_train_7627289.log` / `logs/ab_oas_rod82_7627919.log` / `logs/precompute_ab_oas_7627915.log`.

### Track A v2 — SHM Q0 site-rate prior (`ab_oas_1m_v2_shm`)

v1 CDR mut recall on Rod.82 is **0.150** vs Neutral SHM JC69 **0.688** / CoSiNE **0.658**. v1 already had viral mut-recovery levers (λ_mut=12, entropy, CDR mask×5+force, PSSM, mut-aa-emb) but **Q0 is still ESM-2 MLM** (conservation prior). Fitness tilt β=0 on v1 (default); turning β on would likely **hurt** CDR recovery further.

**v2 recipe (one lever, not a viral-lever pile-on):** mix an SHM site-rate prior into Q0 at train **and** sample:

- Hot sites = IMGT CDR mask ∪ reverse-translated AID **WRCH/DGYW**
- Stay logit −`boost=2.0` on hot sites; +`fwr_stay=0.5` on framework (separate CDR vs FWR rates)
- Destination ratios among mutants stay ESM (mutation head still learns AA identity from OAS pairs)
- Same splits / PLM caches (**7627915**); `SKIP_PRECOMPUTE=1`; `REFRESH_CACHES=0`
- Ckpt `checkpoints/ab_oas_1m_v2_shm` (does **not** overwrite `ab_oas_1m_v1` or pathogen `best.pt`)
- Rod.82 samples `treesbm_ab_oas_v2` (does **not** overwrite `treesbm` / `treesbm_ab` / `treesbm_ab_oas`)
- `--resume`, 72h wall; afterok Rod.82 eval
- Hotspot loss weight 5→8 (same CDR mask, slightly stronger L_mut)

**Not chosen:** S5F/Thrifty dest tables (needs NT precompute + netam); copying viral mrs/λ_mut further (already maxed); ESM fitness β>0 (anti-SHM).

| Step | Job | State | Notes |
|---|---:|---|---|
| train | **7720919** | PENDING (Betty 2026-08-20) | `scripts/slurm_ab_oas_v2_shm_train.sh` → `checkpoints/ab_oas_1m_v2_shm` |
| Rod.82 afterok | **7720920** | dependency on 7720919 | `scripts/slurm_ab_oas_v2_shm_rod82.sh` → `treesbm_ab_oas_v2` |

IDs: [`ab_oas_v2_shm_job_ids.json`](ab_oas_v2_shm_job_ids.json). Metrics stay **—** until afterok eval lands. Do not invent.

---

## Track B — DASM-export TreeSBM retrain + Rodriguez 82 rollout (7539619 / 7539827)

**What was done:** export antibody_benchmark trees excluding the 82 held-out Rodriguez families → train `checkpoints/ab_dasm_v1` → rollout on the frozen 82-family benchmark.

| Step | Job | State | Elapsed | Notes |
|---|---:|---|---|---|
| train | **7539619** | COMPLETED | 01:33:36 | `ab_dasm_v1/best.pt` ep62. **Loss collapsed:** mut/rate/cons=0, `br_target_std=0`, val=0.0038 from epoch 2 |
| rollout+eval | **7539827** | COMPLETED | 02:48:15 | 82 families, config `full_ab_treesbm.yaml` max_seq_len=200 |

**What is being tested:** TreeSBM trained on **collapsed DASM PCP / naive-stripped** trees, then rolled out on the **Rodriguez 82 affinity-maturation families** (same freeze as Track C). This is *not* OAS 1M and *not* CDR-entropy training.

| Model | Root→leaf W1 ↓ | Site freq ρ ↑ | Substitution JS ↓ | Co-mutation ρ ↑ | Leaf diversity W1 ↓ |
|---|---|---|---|---|---|
| treesbm (`ab_dasm_v1`) | 11.080 [10.13,12.09] | 0.144 [0.123,0.164] | 0.508 [0.494,0.520] | 0.043 [0.004,0.083] | 14.269 [13.24,15.43] |

Worse W1/diversity than pathogen-ckpt TreeSBM on the same 82 families (Track C). Expected: collapsed train has no mutational signal.

Holdout: the 82 Rodriguez family IDs in `benchmarks/results/tables/ab_family_ids_82.txt` (sample-igg-SC-* / W-*).

---

## Track C — Rodriguez 82 CoSiNE/DASM/Thrifty benchmark (7517208–12)  pathogen TreeSBM ckpt

**What was done:** frozen 82-family N=20 rollout. TreeSBM used **pathogen** `checkpoints/best.pt` (multi-pathogen flu, job 7062213, ep237 val=16.79) — **not** Ab-trained.

| Model | Job | State | Elapsed |
|---|---:|---|---|
| thrifty | 7517208 | COMPLETED | 00:01:16 |
| dasm_thrifty | 7517209 | COMPLETED | 00:03:50 |
| treesbm (pathogen ckpt) | 7517210 | COMPLETED | 02:36:27 |
| cosine | 7517211 | COMPLETED | 00:58:29 |
| evaluate | 7517212 | COMPLETED | 00:00:38 |

**What is being tested:** CoSiNE / DASM+Thrifty / Thrifty / TreeSBM (OOD pathogen) on **affinity-maturation clonal families** (Rodriguez). Primary metrics = distributional rollout (root→leaf W1, site-freq ρ, subst JS, co-mut ρ, leaf-div W1). TreeSBM is forced to observed topology+BL; CoSiNE is unguided Gillespie.

| Model | Root→leaf W1 ↓ | Site freq ρ ↑ | Substitution JS ↓ | Co-mutation ρ ↑ | Leaf diversity W1 ↓ |
|---|---|---|---|---|---|
| dasm_thrifty | 3.799 | 0.444 | 0.277 | 0.011 | 4.419 |
| cosine | 4.686 | 0.419 | 0.277 | 0.006 | 5.523 |
| thrifty | 4.969 | 0.284 | 0.319 | 0.016 | 8.165 |
| treesbm (pathogen ckpt) | 6.746 | 0.100 | 0.431 | 0.040 | 7.752 |

DASM+Thrifty leads. Pathogen TreeSBM is worst on W1 / site-freq / subst JS (OOD).

Paths: `antibody_benchmark/data/processed/benchmark_trees.jsonl` (82) · summaries `antibody_benchmark/results/summary/` (Track C) vs `antibody_benchmark/results_ab_treesbm/summary/` (Track B).

### Table 5 / Table 6 metrics (Track C + Track B) — 2026-08-15 redo

Paper **Table 5** paste (4 cols; Cov@100 = abs Hamming **e=5**). Full file: [`table6_ab_track_c.md`](table6_ab_track_c.md).

| Method | CDR mut ↑ | SHM err ↓ | Term. div ↓ | Cov@100 ↑ |
|---|---:|---:|---:|---:|
| **Neutral SHM** (JC69 NT, obs. topo) | 0.688 | 0.054 | 17.44 | 0.059 |
| **pLM mutation prior** | — | — | — | — |
| **AR tree-edit** (ARTreeFormer N16 prune + JTT) | 0.000 | 0.117 | 19.14 | 0.081 |
| **CTMC (CoSiNE)** | 0.658 | 0.037 | 3.02 | 0.076 |
| **TreeSBM** pathogen `best.pt` | 0.278 | 0.055 | 5.10 | 0.084 |

- **Neutral ≠ Thrifty.** Thrifty job **7517208** is additive context-SHM only.
- **pLM —:** Rod.82 gen still running on login01 nohup (**bash 4143208 / python 4143211**, CPU, `logs/ab_plm_prior_nohup.log`; samples `antibody_benchmark/results/samples/plm_prior`). ~568/1640 (~29/82 fam) at 2026-08-15 14:12 ET; ETA ~0.5–1h. Do not kill; eval after 82×20.
- **AR:** tested end-to-end on Rod.82; CDR=0 is measured (JTT+adapted pool), not a stub dash.
- **TreeSBM paper row:** pathogen Track C (**7517210** / ckpt **7062213**). OAS **7627919** = additive `treesbm_ab_oas`.
- Additive: dasm_thrifty / thrifty / treesbm_ab / treesbm_ab_oas in `table6_ab_track_c.md`.

---

## Quick mix-up guard

| | Track A OAS 1M | Track B ab_dasm_v1 | Track C Rodriguez 82 |
|---|---|---|---|
| Data | OAS Homo_sapiens 1M → 611 clones | DASM PCP export minus 82 | Frozen 82 families |
| Holdout | 170 clones / 10 donors (SRR17818205-heavy) | the 82 | the 82 |
| Mask | IMGT CDR empirical 27/160 | (DASM export) | n/a |
| Entropy | ON | collapsed | n/a |
| Ckpt | `ab_oas_1m_v1/best.pt` (**7627289**; Rod **7627919** DONE) | `ab_dasm_v1` collapsed | pathogen `best.pt` |
| Jobs | 7585306 OK → … → **7627289** FAIL-with-ckpt / pre **7627915** DONE / Rod **7627919** DONE | 7539619 / 7539827 | 7517208–12 |
| Tested | TreeSBM on OAS trees + CDR mask; Rod.82 Track C metrics pasted | TreeSBM on collapsed DASM, eval on 82 | CoSiNE/DASM vs OOD TreeSBM on 82 |

