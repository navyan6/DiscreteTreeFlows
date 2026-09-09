# S+H26-inspired mut_recovery iterations (TreeSBM kept)

Primary KPI: **COVID `mut_recovery`** (`positional_recovery`: fraction of root≠GT sites where gen==GT).  
Floor: **cons_retention ≥ ~0.98**. Secondary: EVEscape enrichment, Coverage@K / mut F1 (`track_a.py`).

**Constraint:** keep TreeSBM formulation — bridge matching, `log R_θ = log R0 + c_θ`, `RateHeads`. Do **not** replace with antiGen’s standalone mutation classifier.

Wave-1 context (`WAVE1_RESULTS_ANALYSIS.md`): covid_v3_cons at mrs=0.1 → mut≈0.001, cons≈0.97 (over-conserves). Agent G already owns mrs∈{0.3,0.5,1.0} sweeps + v4 (`λ_mut=8`, pos-emb, `L_semi`) eval — coordinate, don’t duplicate those jobs.

antiGen (Specht & Hie 2026) trains on **parent→child phylogenetic mutations**, uses a **joint L×A softmax** (site propensity preserved), and often **gates** with a train PSSM. Map those ideas into TreeSBM rates / losses / sampling below.

---

## Top 5 (do these first)

| # | Idea | Repo change | mut ↔ cons | Effort | Priority | Status |
|---|------|-------------|------------|--------|----------|--------|
| 1 | **Raise eval mutation mass (mrs / steps)** before more arch | `slurm_inf_sweep.sh` / `slurm_eval_covid.sh` @ mrs 0.3–1.0, 12–24h, `qos=mig-max`. Pick best under cons≥0.98. | ↑mut, ↓cons risk at high mrs — stop when cons breaks floor | Low (ops) | **P0** | ✅ **Implemented** (sweep scripts + SKIP_EXISTING; enrichment accepts mrs cleanly) |
| 2 | **Mutation-only / mutation-normalized bridge term** (antiGen loss ∝ NLL of *mutations*, /Z) | `src/bridge/losses.py`: `mut_normalize=count` + `λ_mut` / `λ_cons` + optional `--entropy-weight-alpha-cons`. | Strong ↑mut; cons may drop | Low–med | **P0** | ✅ **Implemented** |
| 3 | **Global site-propensity during sampling** | Inference: `--site-softmax-sample` in `eval_single_tree` / `generate_tree` / enrichment. | ↑mut_recovery via hot sites | Med | **P1** | ✅ **Implemented** |
| 4 | **Static PSSM gate on log rates** | Train+eval: `--pssm-gate` → `log R_eff = w·Z(log R_θ)+(1−w)·Z(log PSSM)`. | ↑mut on variable sites | Med | **P1** | ✅ **Implemented** (default off) |
| 5 | **Per-position sequence features for c_θ** | `--mut-aa-emb` into RateHeads mut head. | ↑AA placement at mut sites | Med | **P1** | ✅ **Implemented** (default off; old ckpts load) |

---

## Flags cheat-sheet

### Train (`scripts/train.py`)

| Flag | Default | Effect |
|------|---------|--------|
| `--lambda-mut` | `5.0` | Scale on `L_mut` inside `L_rate` |
| `--lambda-cons` | `1.0` | Scale on `L_cons` (`L_rate = λ_mut L_mut + λ_cons L_cons`) |
| `--mut-normalize {mean,count}` | `mean` | `count` → Σ(w·kl)/n_mut (stronger hotspot mass) |
| `--entropy-weight-alpha-cons` | same as `--entropy-weight-alpha` | Lower to ease over-conservation |
| `--pssm-gate` | off | Blend Z(log R_θ) with train log-PSSM; learnable per-site `w=σ(γ)` |
| `--pssm-gate-fixed-w` | none | Fixed w∈[0,1] instead of learnable γ |
| `--mut-aa-emb` | off | Concat current-AA emb into c_θ input (new ckpt) |
| `--mut-aa-emb-dim` | `16` | AA emb width |
| `--per-site-pos-emb` / `--lambda-semi` / `--deep-mut-head` | existing | Unchanged |

### Inference (`eval_single_tree` / `generate_tree` / `eval_evescape_enrichment`)

| Flag | Default | Effect |
|------|---------|--------|
| `--mutation-rate-scale` / mrs | `1.0` / `0.3` (enrichment) | Scales mutation fire rate |
| `--site-softmax-sample` | off | Site-propensity categorical → AA\|site |
| `--site-temperature` | `1.0` | Softmax temp on site scores |

Checkpoint `config` + `log_pssm` / `col_entropy` are loaded automatically; old v3/v4 state_dicts still load when new flags are off.

---

## Fast path (covid_v5 + high-mrs)

```bash
# 1) Sync code to Betty, then train v5 (stronger mut + PSSM + AA emb + pos-emb + L_semi)
sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh
# → checkpoints/covid_v5_mutrec  (λ_mut=12, λ_cons=0.5, mut_normalize=count,
#    entropy_alpha_cons=1.0, --pssm-gate --mut-aa-emb --per-site-pos-emb --lambda-semi 0.05)

# Optional ablations:
#   LAMBDA_MUT=16 LAMBDA_CONS=0.25 sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh
#   V4_COMPAT=1 CKPT_DIR=checkpoints/covid_v4_mutrec sbatch --qos=mig-max scripts/slurm_covid_train_mut_recovery.sh

# 2) High-mrs enrichment (skips finished JSONs; fill gaps only)
EVAL_MRS="0.3 0.5 1.0" SKIP_EXISTING=1 sbatch --qos=mig-max scripts/slurm_eval_covid.sh
EVAL_MRS="0.3 0.5 1.0" SKIP_EXISTING=1 sbatch --qos=mig-max scripts/slurm_eval_covid_v4.sh

# 3) Full mrs × n_steps grid (24h) — Agent G may already own v3; check squeue first
sbatch --qos=mig-max scripts/slurm_inf_sweep.sh
# Site-propensity ablation on existing ckpt (no retrain):
SITE_SOFTMAX=1 CKPT=checkpoints/covid_v3_cons/best.pt sbatch --qos=mig-max scripts/slurm_inf_sweep.sh

# 4) After v5 best.pt exists:
CKPT=checkpoints/covid_v5_mutrec/best.pt sbatch --qos=mig-max scripts/slurm_inf_sweep.sh
SITE_SOFTMAX=1 CKPT=checkpoints/covid_v5_mutrec/best.pt sbatch --qos=mig-max scripts/slurm_inf_sweep.sh
```

---

## Full brainstorm (mapped)

### A. Training objective / data (antiGen: mut NLL on phylo edges)

| Idea | What to change | mut ↔ cons | Effort | Priority |
|------|----------------|------------|--------|----------|
| A1. Mutation-normalized L_mut (Top-5 #2) | `losses.py` + train flags | ↑mut / ↓cons risk | L–M | P0 ✅ |
| A2. Downweight L_cons / entropy-cons | `--lambda-cons`, `--entropy-weight-alpha-cons` | ↑mut / ↓cons | L | P0 ✅ |
| A3. Explicit parent→child mut auxiliary | From each bridge edge, add CE/NLL on observed substitutions (antiGen-style) as `L_edge_mut` next to bridge KL; λ_edge small | ↑mut specificity | M | P2 |
| A4. De-novo / rare-mut upweight | Weight mut positions by inverse train mutation frequency (or 1 for never-seen AA) inside L_mut | ↑hard muts; may hurt cons slightly | M | P2 |
| A5. Temporal COVID train split | Train on `data/covid_temporal` (see `SPLITS.md`) so objective matches forecasting; eval mut_recovery on post-cutoff leaves | Fairer ↑mut for Table 3/5 narrative | M (data+train) | P2 (after geo mut moves) |

### B. Architecture (keep log R0 + c_θ)

| Idea | What to change | mut ↔ cons | Effort | Priority |
|------|----------------|------------|--------|----------|
| B1. AA + pos into mut head (Top-5 #5) | `RateHeads` `--mut-aa-emb` | ↑mut | M | P1 ✅ |
| B2. Deep mut head | Already: `--deep-mut-head` / `DEEP_MUT=1` → new ckpt dir (`covid_v5_deepmut`) | Mild ↑mut capacity | L | P2 if v4 flat |
| B3. Mask identity AA in rate softmax | When converting rates→probs for sampling, set stay-AA logit −∞ (antiGen zeros non-mutations) so mass goes to true substitutions | ↑mut events; cons ↓ unless site fire rate controlled | L | P1 (partially via site-softmax mut dists) |
| B4. Finetune small linear on log_R0 toward mut NLL | Freeze GraphTF; train only mut-head (+ optional R0 scale) on mut positions | ↑mut with less overfitting risk | M | P2 |
| B5. Larger ESM R0 (playbook step 5) | Swap ESM size in precompute / train | Uncertain; Table 7 direction | H | P3 |

### C. Inference / hybrid (antiGen: PSSM gate, recall@q)

| Idea | What to change | mut ↔ cons | Effort | Priority |
|------|----------------|------------|--------|----------|
| C1. mrs / n_steps / branch_rate_scale grid (Top-5 #1) | Existing sweep scripts; also tune `branch_rate_scale` | ↑mut | L | P0 ✅ |
| C2. Site-propensity sampling (Top-5 #3) | `--site-softmax-sample` | ↑mut | M | P1 ✅ |
| C3. PSSM / EVEscape static gate (Top-5 #4) | `--pssm-gate` | ↑mut + enrichment | M | P1 ✅ |
| C4. Report next-mutation recall / Coverage@K | Align paper metric with antiGen recall@q via `track_a.py` + optional ranking of `log R_θ` over L×20 vs GT mut set | Diagnostic (doesn’t change gen) | L | P1 (metric) |
| C5. Multi-sample best-of-M leaf | Generate M leaves, pick max mut_recovery under cons floor | ↑reported mut (honest if captioned) | L | P2 |

### D. Eval protocol hygiene (don’t confuse KPI)

| Idea | What to change | mut ↔ cons | Effort | Priority |
|------|----------------|------------|--------|----------|
| D1. Finish v3/v4 enrichment JSONs | Agent G jobs; do not invent numbers | — | Ops | P0 |
| D2. Debug H1N1 leaf-holdout identity 0.08 | `eval_leaf_holdout.py` / matching before trusting Table 3 | — | M | P1 (Table 3, not COVID mut) |
| D3. Clade-holdout COVID | `data/covid_cladeholdout` + same mut_recovery | Secondary COVID KPI | M | P2 |

---

## Suggested decision tree (post–Agent G sweeps)

```
mrs grid on v3/v4 done?
  ├─ mut rises ≥ BD/PLM and cons≥0.98 → STOP arch; fill Table 5
  ├─ mut rises but cons < 0.98 → lower mrs or add PSSM gate / cons λ
  └─ mut still ~0 at mrs≥0.5
        → A1/A2 (λ_mut↑, λ_cons↓) retrain v5   ← covid_v5_mutrec fast path
        → then B3 + C2 (sampling) without retrain
        → then B1 / C3 if still flat            ← also on by default in v5 train
```

## Explicit non-goals

- Replacing bridge matching with antiGen’s standalone Transformer classifier.
- Claiming antiGen numbers or inventing Table 2/3/5 cells.
- Conflicting with Agent G’s in-flight `eval_covid` / `inf_sweep` / `covid_v4_mutrec` jobs — extend playbook here; submit only baselines / new ablations not already queued.
