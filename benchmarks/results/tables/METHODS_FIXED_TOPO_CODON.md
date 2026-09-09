# Methods update — fixed topology, codon Q0, NT stall fire

Use this for the paper methods (viral sequence model). **No new numbers here** until Betty jobs finish. Paper checkpoints were not overwritten.

## Motivation

`aa|hit` is destination identity given a substitution. TreeSBM’s default generator is independent-site AA with an ESM residual. Neutral CTMC + BD and PhyloFlow look strong on coverage / RF in part because (i) COVID Cov@ε=2 saturates, (ii) Neutral uses a real \(Q\), (iii) PhyloFlow trains topology. This wave isolates **destination** from **topology**.

## 1. Fixed-topology sequence evaluation (cheapest; no train)

**Protocol.** For each held-out test tree, keep the **observed Newick and branch lengths**. Paint sequences from the true root down the tree. Score leaves **by node id** (not nearest-neighbor matching).

Primary KPIs (unchanged definitions): `mut_recovery`, `site_recall`, `aa_acc_given_hit`, `cons_retention`.

Painters (same independent-site Bernoulli on observed BL; only \(Q\) / residual changes):

| Painter | Substitution process |
|---|---|
| Neutral | Flat 20×20 AA CTMC \(Q\) (same generator as NeutralBD; **no free birth–death**) |
| Codon-GY94 | Goldman–Yang codon \(Q\) (κ=2, ω=0.25), uniform mixture over sense codons, collapsed to 20 AA |
| TreeSBM | Frozen paper ckpt RateHeads residual \(c_θ\) on its trained R0 (ESM-2-8M), parent-only tree context per edge |

This is the viral analogue of the antibody fair protocol (forced observed topology). It answers whether Neutral/PhyloFlow “win because they learn topology.” Topology is identical across rows.

**Script:** `scripts/eval_fixed_topo_seq.py`. **Jobs:** `scripts/slurm_eval_fixed_topo.sh` (COVID / H1 / H3 / HIV geo paper ckpts, read-only).

## 2. Codon Q0 + learned residual (train; new dir)

**Model.** Same TreeSBM architecture as COVID `covid_v5_mutrec` (site pos-emb, empirical entropy weighting). Replace frozen R0 with GY94→AA (`--r0-backend codon_gy94`). Residual \(c_θ\) is still learned. Fitness tilt **β=0**.

Caches: `group_*_ref_rates_codon_gy94.pt` (does **not** overwrite untagged ESM `group_*_ref_rates.pt`).

Checkpoint: `checkpoints/covid_v8_codon_q0/` (new). Locked dirs (`covid_v5_mutrec`, flu/HIV paper, OAS v1–v3, mutlin) are refused.

After train, the same fixed-topology eval is run on `covid_v8_codon_q0/best.pt`.

## 3. NT stall flags on fire (eval-only; no dest change)

Homopolymer (≥4), G-quadruplex, palindrome, and tandem-repeat flags are computed on a **preferred-codon reverse-translate** of the current AA (not authentic CDS). On flagged sites, **fire** \(1-p_{\mathrm{stay}}\) is multiplied by \(1+\alpha\) (default α=1); destination ratios among non-parent AAs are unchanged.

Applied at paint time as extra rows `*_ntfire`. Not written into paper checkpoints.

## What this is not

- Not an indel / template-switch model.
- Not free topology (BD or PhyloFlow).
- Antibody Recipe C (cosine RateHeads) is unchanged.

## Jobs

IDs: `fixed_topo_codon_job_ids.json`. Free-gen enrich: `covid_v8_eval_job_ids.json`.

## Results — `covid_v8_codon_q0` fixed-topo (pulled 2026-08-26)

Train **7828383** DONE; fixed-topo eval **7828384** DONE (`covid_v8_codon_q0_fixed_topo.json`).  
Ckpt: `checkpoints/covid_v8_codon_q0/best.pt` (eval log: epoch **30**, val **2.1505**). Protocol: `data/covid/test`, mrs=1.0, N=20, paint on observed Newick+BL.

| Painter | aa\|hit ↑ | mut_recovery | site_recall | cons_retention ↑ |
|---|---:|---:|---:|---:|
| Neutral | 0.002 | ~0 | ~0 | 0.999 |
| Codon-GY94 | 0.079 | ~0 | 0.0078 | 0.999 |
| **TreeSBM v8** (GY94 R0 + residual) | **0.089** | ~0 | ~0 | 0.999 |
| TreeSBM v8 + NT fire | **0.139** | ~0 | ~0 | 0.998 |
| Codon-GY94 + NT fire | 0.086 | ~0 | ~0 | 0.999 |

**vs paper `covid_v5_mutrec` same protocol** (`covid_fixed_topo_mrs1.0.json`): v5 TreeSBM aa\|hit **0.113** (ntfire **0.320**). **v8 does not beat v5** on fixed-topo aa\|hit; GY94 R0 alone (0.079) is close to v8 residual (0.089). Mut/site recall stay near zero under short observed BLs — dest bake-off only.

Free-gen enrichment mrs=0.5 (**7884601** DONE): `eval_enrichment_covid_v8_codon_q0_mrs0.5.json` → aa\|hit **0.087**, mut_recovery **0.032**, site_recall **0.283**, cons **0.729**. Paper `covid_v5_mutrec` free-gen aa\|hit ≈ **0.375** — **codon GY94 R0 collapses destination identity** under free gen, not just fixed-topo.

### Why codon / “chemical” sequence features hurt aa\|hit

`aa|hit` = P(dest AA = GT | both left parent). It only cares about **which of 19 mutants** you pick, not how often you leave.

1. **GY94→AA R0 is a wrong dest prior for Spike.** Codon GY94 (κ, ω) encodes purifying selection + transition bias in codon space, then we collapse to AA. Observed SARS-CoV-2 Spike evolution is immune-escape / epistatic / site-heterogeneous — not a single ω CTMC. ESM-2 MLM residual (v5) was already a better dest prior than GY94 (fixed-topo 0.113 vs 0.079; free-gen ≈0.375 vs **0.087**).
2. **Fire knobs ≠ dest knobs.** NT stall fire only scales \(1-p_{\mathrm{stay}}\); dest ratios among non-parent AAs are unchanged. Fixed-topo COVID ntfire *raises* aa\|hit (v5 0.11→0.32) by changing *which sites fire* under short BLs — that is sampling noise / site selection, not learning identity. It does not fix a bad R0.
3. **Stall/chem biology is the wrong channel for this metric.** Homopolymers/G4/hairpins drive indels & template switches; the sampler is AA substitution only. Stratified audit: stall windows are *not* where aa\|hit collapses; chem\|hit − aa\|hit ≈ 0.09 means we often land BLOSUM-neighbors, not the GT residue — soft chemical similarity without exact identity.
4. **Training against GY94 R0 dilutes the residual.** v8 residual barely beats bare GY94 (0.089 vs 0.079 fixed-topo) — the network did not unlearn the codon prior enough to recover ESM-level dest.
