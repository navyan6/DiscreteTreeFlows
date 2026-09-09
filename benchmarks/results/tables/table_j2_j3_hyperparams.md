# Appendix J.2 / J.3 — Hyperparameters & sampling settings

**Status 2026-08-17:** LaTeX filled for `covid_v5_mutrec` (`tab:app_hyperparameters`, `tab:app_sampling_settings`). Paste: [`table_j2_j3_hyperparams.tex`](table_j2_j3_hyperparams.tex).

**Ckpt:** `checkpoints/covid_v5_mutrec/best.pt` — best **epoch 14** / val=0.9931 (train job **7350224**, timed out epoch 42; configured max 100, patience 50).

Sources: Betty ckpt + `logs/covid_mutrec_7350224.log`, `scripts/slurm_covid_train_mut_recovery.sh`, `scripts/train.py`, `src/networks.py`, `scripts/slurm_table5_coverage.sh`, `benchmarks/generation.py`, `benchmarks/results/tables/TABLE_METADATA.md`, `results/wave5_pull/ckpt_meta/covid_v5_mutrec_meta.json`.

---

## J.2 Model hyperparameters — paste-ready

| Hyperparameter | Value | Source |
|---|---|---|
| Tree encoder layers | **4** | `scripts/train.py` (`TreeEncoder n_layers=4`) |
| Hidden dimension | **128** | `scripts/train.py` (`d_model=128`) |
| Attention heads | **8** | `scripts/train.py` (`n_heads=8`) |
| pLM reference | **ESM-2 8M** (`facebook/esm2_t6_8M_UR50D`, d_plm=320) | `src/r0_backends.py` / train |
| Mutation head depth | **2** (Linear→ReLU→Linear; `deep_mut_head=false`) | `src/networks.py` + ckpt meta |
| Branching head depth | **2** (d_model→64→1) | `src/networks.py` |
| Branch-length head depth | **2** (d_model→64→1) | `src/networks.py` |
| Dropout | **0.1** | `scripts/train.py` |
| Optimizer | **AdamW** (weight_decay=1e-4; cosine `T_max=epochs`, `eta_min=1e-6`) | `scripts/train.py` |
| Learning rate | **1e-4** | `scripts/train.py` / slurm |
| Batch size | **1** tree | `scripts/train.py` DataLoader |
| Training epochs | **14** (best ckpt; configured max **100**, patience **50**) | ckpt + job **7350224** log |
| Max seq len (COVID) | **1280** | `slurm_covid_train_mut_recovery.sh` |
| Bridge constant `c` | **1.0** | `--bridge-c` |
| λ_mut / λ_cons / λ_br / λ_top | **12.0** / **0.5** / **0.1** / **0.1** | slurm + meta |
| λ_semi | **0.05** | slurm |
| Site entropy weighting | ON; α=**3.0**, floor=**1.0**, source=`empirical` | meta |
| Per-site pos emb / mut AA emb / PSSM gate | ON; d_aa=**16** | meta |

---

## J.3 Sampling settings — paste-ready

Locked to **Table 5 Coverage@K=100** protocol (ckpt `covid_v5_mutrec`, job **7575015**). Enrichment companion (mut recall / EVEscape) uses **n_steps=100**, **mrs=0.5**, **max_trees=20**, **seed=42** — not all appear in the paper J.3 table rows.

| Setting | Value | Source |
|---|---|---|
| Evolutionary horizon (`n_steps`) | **50** | `slurm_table5_coverage.sh` / `coverage_curves.py` |
| Maximum nodes (`max_leaves`) | **400** | `benchmarks/generation.py` default; E.3 prod row |
| Maximum branching factor | **2** (binary bifurcation: 0 or 2 children) | `eval_single_tree.py` / `generate_tree.py` |
| Mutation temperature (`site_temperature`) | **1.0** | default; E.3 sweep flat (**7639276–77**) |
| Branching temperature | **1.0** (learned $\lambda$; scale **6.0**) | no branching-temp knob; `branch_rate_scale=6.0` at inference (`generate_tree.py`) |
| Stop threshold | **none** (fixed `n_steps`) | `stop_prob` head trained but **not applied** at inference |
| K (Coverage@100) | **100** | `slurm_table5_coverage.sh` (`K_MAX=100`) |
| Seeds | **0** (coverage roots); enrichment uses **42** | `coverage_curves.py` / `eval_evescape_enrichment.py` |
| `branch_rate_scale` | **6.0** | `generate_tree.py` / `TreeSBMMethod` |
| `mutation_rate_scale` (enrichment) | **0.5** | TABLE_METADATA / T5/T7/T8 |
| N (held-out topology) | **16** | coverage / Table 5 |
| `(N,H)` adapter | `rate_per_H=1.2`, `cushion=1.6`, `max_retries=4` | `benchmarks/methods/treesbm.py` |
| Enrichment `n_steps` | **100** (when scoring mut/cons/EVEscape) | `slurm_eval_covid.sh` / TABLE_METADATA |
| Enrichment `max_leaves` | **300** (script default) | `eval_evescape_enrichment.py` |
| Coverage adapter cap @ N=16 | **28** (= ⌊1.6·16⌋+2) unless `--max-leaves` override | `TreeSBMMethod` |

**Note:** E.3 temp ∈ {0.5,1.0,1.5} and max-leaves ∈ {100,400,800} did not move Cov@e2 on locked Brazil roots — see [`table_e3_sampling.md`](table_e3_sampling.md).
