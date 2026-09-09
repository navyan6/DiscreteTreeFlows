# Antibody affinity-maturation forward-generation benchmark

Compares **TreeSBM, CoSiNE, DASM+Thrifty, and Thrifty** on **root-conditioned recursive rollout** over held-out clonal trees (fixed observed topology + branch lengths).

See [`DATA_AUDIT.md`](DATA_AUDIT.md) before changing data sources.

## Quick smoke (null models)

```bash
# from repo root
python antibody_benchmark/scripts/build_dataset.py --config antibody_benchmark/configs/smoke.yaml
python antibody_benchmark/scripts/run_rollouts.py --config antibody_benchmark/configs/smoke.yaml
python antibody_benchmark/scripts/evaluate.py --config antibody_benchmark/configs/smoke.yaml
```

## Full models (Betty)

```bash
bash antibody_benchmark/scripts/download_data.sh
pip install -e antibody_benchmark/data/raw/repos/netam
# CoSiNE: cd antibody_benchmark/data/raw/repos/cosine && uv sync

# Edit configs/default.yaml — set treesbm.checkpoint; keep cosine.guided=false
python antibody_benchmark/scripts/build_dataset.py --config antibody_benchmark/configs/default.yaml
python antibody_benchmark/scripts/run_rollouts.py --config antibody_benchmark/configs/default.yaml
python antibody_benchmark/scripts/evaluate.py --config antibody_benchmark/configs/default.yaml
```

## Fairness (non-negotiable)

1. Same held-out families / root / topology / BLs / N rollouts  
2. Generated parents feed children (no teacher-forcing of true internals)  
3. Family-level train/test isolation  
4. CoSiNE: unguided Gillespie only  
5. TreeSBM: forced observed topology+BL for this experiment  

Neutral SHM (Thrifty Dryad OOF/synonymous) is a **separate** optional benchmark.
