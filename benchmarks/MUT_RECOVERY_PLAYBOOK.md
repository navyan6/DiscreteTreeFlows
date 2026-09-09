# Mut-recovery architecture iteration playbook (Agent G)

Primary KPI: **mut_recovery**. Floor: **cons_retention ≥ ~0.98**.
Secondary: Coverage@K / mutation-set F1 via `benchmarks/track_a.py`.

## Iteration order

1. **Pos emb + entropy (already in covid_v3_cons path)**  
   Flags: `--per-site-pos-emb --use-site-entropy --use-entropy-loss-weighting --use-entropy-cons-weighting --entropy-source empirical`  
   Train: `scripts/slurm_covid_train_cons.sh` (v3) or `scripts/slurm_covid_train_mut_recovery.sh` (v4 = +pos emb + L_semi + λ_mut≥8).

2. **Stronger mutation head / λ_mut**  
   Env: `LAMBDA_MUT=8.0` (default in mut_recovery SLURM).  
   Deeper MLP: `DEEP_MUT=1 CKPT_DIR=checkpoints/covid_v5_deepmut` → `--deep-mut-head` (mut_in→128→64→20). **Do not** resume into v4 with this flag.

3. **Enable L_semi**  
   `--lambda-semi 0.05` (default in mut_recovery SLURM). Ablate `{0, 0.01, 0.05, 0.1}`.

4. **Inference sweep**  
   `sbatch --qos=mig-max scripts/slurm_inf_sweep.sh`  
   Defaults: `mutation_rate_scale ∈ {0.3,0.5,1.0}` × `n_steps ∈ {50,100,150}` (skips finished tags; 24h wall).  
   Full script: `scripts/slurm_inference_sweep.sh`. COVID smoke used 0.3 / 100.

5. **Optional larger ESM-2 for R0** (Table 7 direction) — only if still behind PLMPrior.

## Stop criterion

TreeSBM `mut_recovery` > NeutralBD / PLMPrior / adapted topology methods on the **same** eval protocol (`eval_evescape_enrichment.py` or `eval_leaf_holdout.py` / `track_a.py`).

## Eval commands (Betty, qos=mig-max)

```bash
# COVID mut/cons + EVEscape — v3 @ mrs 0.3/0.5 (12h)
sbatch --qos=mig-max scripts/slurm_eval_covid.sh

# High-mrs inference sweep (24h; skips existing)
sbatch --qos=mig-max scripts/slurm_inf_sweep.sh

# covid_v4_mutrec @ mrs 0.3/0.5 (when best.pt exists)
sbatch --qos=mig-max scripts/slurm_eval_covid_v4.sh
# or one-liner on Betty:
bash scripts/run_eval_covid_v4_ready.sh

# H1N1 leaf-holdout
sbatch --qos=mig-max scripts/slurm_h1n1_leafholdout_eval.sh

# Track A (Coverage@K)
python benchmarks/track_a.py --checkpoint checkpoints/covid_v4_mutrec/best.pt \
    --data data/covid/test --K 10 --n-groups 20 \
    --out benchmarks/results/track_a_covid_v4.json
```

## Submitted Agent G jobs (2026-08-01)

| JobID | Script | Notes |
|------:|--------|-------|
| 7349962 | `slurm_eval_covid.sh` | v3 mrs 0.3/0.5 |
| 7349963 | `slurm_inf_sweep.sh` | FAILED (path bug); superseded |
| **7350021** | `slurm_inf_sweep.sh` | v3 mrs 0.3/0.5/1.0 × steps (resubmit) |
| 7350018 | `slurm_eval_covid_v4.sh` | v4 best.pt @ mrs 0.3/0.5 |
| 7323331 | (prior) `covid_mutrec` | still training; ~9.7h wall left at submit time |

## Decision tree after jobs land

1. If **v3 mrs≥0.3** lifts mut with cons≥~0.98 → lock sampling; use for Table 5; compare v4.  
2. Else if **v4** lifts mut → prefer v4; optional `CKPT=.../covid_v4_mutrec/best.pt sbatch --qos=mig-max scripts/slurm_inf_sweep.sh`.  
3. Else → `LAMBDA_MUT=12` and/or `DEEP_MUT=1` into a **new** ckpt dir; debug leafholdout identity=0.08 before Table 3.
