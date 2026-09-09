# Tree-generation table — runbook

Fills "Evolutionary tree generation on held-out roots" (8 methods × 6 metrics),
two tracks (simulated + empirical). Native rows run in the treesbm env; adapted
rows need the external repos (see `EXTERNAL.md`); PhyloGFN/BHV-flow are excluded
(see `BLOCKERS.md`).

## 0. deps (treesbm env, on the cluster)
```bash
pip install dendropy pyvolve tqdist
```

## 1. fit BD/substitution params on TRAIN only
```bash
python benchmarks/fit_params.py --train-data data/h3n2/train --out benchmarks/results/params.json
```

## 2. leakage check (train ∩ test node ids = ∅)
```bash
python benchmarks/heldout/build_examples.py --check-leakage \
    --train-data data/h3n2/train --test-data data/h3n2/test
```

## 3. dev run (native 4, N=16, tiny) — sanity end-to-end
```bash
python benchmarks/run_table.py --test-data data/h3n2/test \
    --params benchmarks/results/params.json \
    --checkpoint checkpoints/h3n2_temporal/best.pt \
    --N 16 --K 10 --M 20 --max-roots 5 --out benchmarks/results/results_dev.csv
python benchmarks/make_table.py --results benchmarks/results/results_dev.csv \
    --out benchmarks/results/tables_dev
```

## 4. full run (native 4)
```bash
python benchmarks/run_table.py --test-data data/h3n2/test \
    --params benchmarks/results/params.json \
    --checkpoint checkpoints/h3n2_temporal/best.pt \
    --N 16 32 64 --K 100 --M 50 --max-roots 100 --out benchmarks/results/results.csv
python benchmarks/make_table.py --results benchmarks/results/results.csv --ci boot
```

## 5. adapted rows (ARTreeFormer, PhyloVAE)
Produce topology pools per `EXTERNAL.md`, then extend `build_methods` in
`run_table.py` with two `TopologyPriorMethod` instances (pool + shared adapters).

## unit tests (pure-Python; run anywhere)
```bash
for t in heldout metrics_trees metrics_seq_dist metrics_matching \
         metrics_matched validity_table events_taskdata; do
  python benchmarks/tests/test_$t.py; done
```

## Structure
```
benchmarks/
  heldout/build_examples.py   held-out roots + induced subtrees (+ collapse, H, leakage)
  fit_params.py               BD/substitution/AA-freq fit on train (Yule)
  methods/                    base, bd_topology, bd_methods, plm_prior, treesbm, topology_prior
  adapters/                   branch_length, sequence (pyvolve / ESM)
  sim/reference.py            simulated reference distributions
  metrics/                    trees, distributions(+JS), sequences, matching, matched, branch_lengths, events
  validity.py                 per-sample validity checks
  run_table.py                orchestrator -> long-format results.csv
  make_table.py               CSV + LaTeX (mean ± SE / bootstrap CI)
  EXTERNAL.md / BLOCKERS.md    external-repo procedure + honest exclusions
```
