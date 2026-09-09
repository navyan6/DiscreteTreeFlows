# Antibody benchmark — freeze + smoke status

**Date:** 2026-08-10 (Betty Priority-6 smokes complete)

## Frozen artifacts (done)

| Artifact | Path | Status |
|---|---|---|
| Trees JSONL | `antibody_benchmark/data/processed/benchmark_trees.jsonl` | **82** families |
| Family IDs | `antibody_benchmark/data/processed/benchmark_family_ids.txt` | 82 lines |
| Attrition CSV | `antibody_benchmark/data/processed/data_audit.csv` | 7769 rows |
| Filter summary | `antibody_benchmark/data/processed/FILTER_SUMMARY.md` | written |
| Held-out summary | `antibody_benchmark/data/processed/heldout_summary.json` | written |

Source PCP md5 `05719ad05f8e1e8b3bf99a8daf9b76ff` (Betty canonical DASM; not redownloaded).

## Multi-root / naive policy

- Zenodo v3 PCPs are upstream `no-naive` → multi-root forests after germline edge removal.
- Naive retention implemented; Rodriguez has 0 `parent_is_naive` → all 82 use `root_source=topological`.
- Attrition: 7769 → 7191 after root → **82** eligible.

## Unit tests

`pytest antibody_benchmark/tests/test_reconstruct_and_rollout.py` → **9 passed**.

## Model smokes — Betty (`treesbm` env, family `sample-igg-SC-24::1117`)

| Model | Result | Notes |
|---|---|---|
| Thrifty | **PASS** | 7 nodes; 2026-08-10T14:17Z |
| DASM+Thrifty | **PASS** | 7 nodes; 2026-08-10T14:18Z |
| TreeSBM | **PASS** | 7 nodes; ckpt `~/DiscreteTreeFlows/checkpoints/best.pt` |
| CoSiNE | **PASS** | 7 nodes; 2026-08-10T19:08Z on B200 MIG (`flash-attn` 2.8.3.post1) |

Reports: `antibody_benchmark/results/smoke_reports/{thrifty,dasm_thrifty,cosine,treesbm}.json`.

### Local Mac (unchanged)

Thrifty/DASM **PASS**; CoSiNE **FAIL** (no CUDA flash-attn); TreeSBM **FAIL** (no local `checkpoints/best.pt`).

## Checkpoint paths (Betty)

| Artifact | Path | Verified |
|---|---|---|
| CoSiNE | `/vast/projects/pranam/lab/nnori/antibody_benchmark/cosine_ckpts/cosine_dasm.ckpt` | size 1789059027; md5 `42b14ca6b41e0689f46eb17b0b52eb65`. Symlinked into repo checkpoints dir. |
| TreeSBM | `/vast/home/n/nnori/DiscreteTreeFlows/checkpoints/best.pt` | real Betty artifact (Jul 12; keys `node_enc`/`tree_enc`/`rate_heads`) |

## Betty ops notes

- SSH via ControlMaster mux (BatchMode OK while mux live).
- Synced lean code + ckpt via rsync/scp (no git push).
- Built `flash-attn==2.8.3.post1` with `FLASH_ATTN_CUDA_ARCHS=100` in `treesbm` env.
- CoSiNE GPU smoke: `/vast/projects/pranam/lab/nnori/antibody_benchmark/cosine_smoke.sbatch` (job 7512532).

## Remaining blockers

1. Dryad neutral SHM — **DEFERRED** (optional secondary).
2. Local Mac still cannot run CoSiNE/TreeSBM smokes (expected).

## Full benchmark submit (Betty)

**Submitted:** 2026-08-11T00:13:43Z  
**Config:** `antibody_benchmark/configs/full.yaml` — **82 trees × N=20**  
**Job IDs file:** `/vast/projects/pranam/lab/nnori/antibody_benchmark/full_submit_ids.json` (local copy: `results/smoke_reports/full_submit_ids.json`)

| Model | JobID | Partition | Timelimit | State at submit |
|---|---|---|---|---|
| Thrifty | **7517208** | genoa-std-mem | 24h | RUNNING |
| DASM+Thrifty | **7517209** | genoa-std-mem | 24h | RUNNING |
| TreeSBM | **7517210** | genoa-std-mem | 48h | RUNNING |
| CoSiNE | **7517211** | b200-mig45 (qos=mig-max) | 48h | PENDING (ReqNodeNotAvail) |
| evaluate (afterok all 4) | **7517212** | genoa-std-mem | 4h | PENDING (Dependency) |

**Outputs**
- Samples: `~/DiscreteTreeFlows/antibody_benchmark/results/samples/<model>/<family>/rollout_*.json`
- Summary (post-eval): `~/DiscreteTreeFlows/antibody_benchmark/results/summary/`
- Logs: `/vast/projects/pranam/lab/nnori/antibody_benchmark/logs/ab_full_<jobname>_<jid>.{out,err}`

**Submit fixes applied (PARCC CLI filter):** partition `genoa-std-mem` (not `genoa-std`); CPU jobs `--mem-per-cpu=5632M`; `module load slurm` on login.

**Still deferred:** Dryad neutral SHM.
