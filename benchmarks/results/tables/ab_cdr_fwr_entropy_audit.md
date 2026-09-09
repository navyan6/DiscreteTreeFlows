# Antibody CDR vs FWR entropy / mask audit

SHM should elevate mutability in CDR (and AID motifs) vs FWR. If TreeSBM Q0 = ESM MLM, stay mass concentrates on conserved FWR — opposite of Neutral/Thrifty/CoSiNE.

| set | n_seqs | L | n_CDR | n_FWR | H_CDR | H_FWR | H_CDR/H_FWR | mask_frac |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OAS_train_tips | 1873 | 160 | 27 | 133 | 2.2453 | 1.3804 | 1.627 | 0.169 |
| Rod82_leaves | 964 | 160 | 27 | 133 | 2.0540 | 1.0149 | 2.024 | 0.169 |

CDR mask: `results/ab_cdr_mask/mut_hotspot_mask_imgt_cdr_1m.pt`

## Interpretation

- If H_CDR/H_FWR ≫ 1 on OAS/Rod, empirical SHM load is CDR-biased.
- v1 TreeSBM CDR recall 0.15 vs JC69 ~0.69 is a Q0 mismatch, not missing λ_mut.
- Recipe A (`--shm-site-boost` / `--shm-fwr-stay` / `--shm-use-aid`) lowers stay on CDR∪AID.

Thrifty/DASM available as `DASMThriftyModel` in `antibody_benchmark/models/dasm.py`: Thrifty μ × AA selection targets NT SHM hotspots; TreeSBM v1 Q0 is ESM conservation (anti-correlated with CDR entropy). Audit H_CDR/H_FWR ≈ **1.63** (OAS) / **2.02** (Rod.82) confirms CDR-biased diversity under the IMGT mask.

