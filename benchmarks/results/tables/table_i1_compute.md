# Appendix I.1 — Computational cost

**Status 2026-08-17:** LaTeX paste in [`table_i1_compute.tex`](table_i1_compute.tex). **GPU peak mem** still **—** (no `nvidia-smi` / sacct MaxRSS artifact).

**Paper columns:** Method, Params, Train time, GPU mem, Time/root, Notes (`tab:app_compute`)

**Protocol (inference):** COVID Brazil Table 5 — N=16, K\_max=100 pool gen, 5 roots (groups 4/5/7/9/10); ckpt `covid_v5_mutrec` @ ep14. Source: `coverage_curves_covid_N16_table5_eabs.csv` `runtime_gen_sec`.

---

## Paste-ready (paper table)

| Method | Params | Train time | GPU mem | Time/root | Notes / artifact |
|---|---|---|---|---|---|
| pLM mutation prior | ESM-2 **8M** (frozen cache) | **---** (no TreeSBM train) | -- | **3 s** | `runtime_gen_sec`=15.2 / 5 roots |
| Autoregressive tree-edit (ARTreeFormer-adapted) | **--** · Blocker: external ckpt not inventoried | **--** · Blocker: AR train wall not in local logs | -- | **31 s** | `runtime_gen_sec`=153.0 / 5 roots |
| TreeSBM (`covid_v5_mutrec`, L=1280) | **11M** total‡ | **8.0 h** wall to best @ **ep14** (job **7350224**) | -- | **71 min** | `runtime_gen_sec`=21421 / 5 roots |

‡ **11M** = **~3M trainable** (NodeEncoder + TreeEncoder + RateHeads; `best.pt` 11.6 MB ≈ 2.9M float32) + **8M frozen** ESM-2 PLM cache (same backbone as pLM row; not trained in TreeSBM step). Paper **Params** column reports total stack at inference.

### Param reconciliation (11M vs ~3M)

| Component | Count | Source |
|---|---:|---|
| ESM-2 `esm2_t6_8M_UR50D` | ~8M | frozen; cached `*_plm.pt` |
| TreeSBM trainable | ~2.9M | `checkpoints/covid_v5_mutrec/best.pt` byte size / 4 |
| **Total (TreeSBM row)** | **~11M** | sum |

---

## Train wall (`covid_v5_mutrec`)

| Field | Value |
|---|---|
| Slurm job | **7350224** (`b200-mig45`, `--time=24:00:00`) |
| Start | 2026-08-01 21:22:45 EDT |
| `best.pt` @ ep14 | 2026-08-02 05:25:12 EDT |
| **Wall to best** | **8 h 2 min** → **8.0 h** in table |
| Job end | Cancelled @ 24 h time limit (continued to ep42; ep14 retained) |
| Log | `logs/covid_mutrec_7350224.log` (Betty) |

---

## Inference timing (COVID Table 5)

| Method | `runtime_gen_sec` | / 5 roots | Time/root (table) |
|---|---:|---:|---|
| plm\_prior | 15.2 | 3.0 s | 3s |
| artreeformer\_adapted | 153.0 | 30.6 s | 31s |
| treesbm | 21421.5 | 4284 s ≈ 71.4 min | 71m |

Pool wall for 5 roots × K≤100 trees (scoring cheaper after pool gen).

---

## Extended rows (not in paper `tab:app_compute`)

| Method | Params | Train time | GPU mem | Time/root | Time/100 trees | Notes |
|---|---|---|---|---|---|---|
| Neutral BD | substitution matrix only | --- | -- | ~38 s (COVID) | ~188 s / 5 roots | same CSV |
| BHV / PhylaFlow | -- | -- | -- | -- | -- | See `EXTERNAL.md` |
| TreeSBM (H3N2 HA L=566) | same arch | `h3n2_v3_lit_hotspot` | -- | ~52 min | ~4.35 h / 5 roots | `runtime_gen_sec`=15668 |
| TreeSBM Ab OAS 1M | same arch, L=160 | **05:15:26** @ ep122 (job **7627289**) | -- | -- | -- | `antibody_tracks.md` |

### C.3 generation job walls (Slurm Elapsed; multi-method)

| Virus | Job | Elapsed |
|---|---:|---|
| covid | 7627191 | 05:57:38 |
| h1n1 | 7627192 | 02:54:29 |
| h3n2 | 7627193 | 04:26:08 |
| hiv_geo | 7627194 | 01:30:07 |
| hiv_temporal | 7627195 | 02:30:36 |

---

## Blockers

| Cell | Blocker |
|---|---|
| GPU mem (all) | No MaxRSS / `nvidia-smi` artifact; sacct unavailable on Betty login path |
| AR Params + train time | External ARTreeFormer baseline; train wall not in repo logs |
| Exact trainable `numel` | Optional torch recount (~2.9M from ckpt size is sufficient for paper) |
