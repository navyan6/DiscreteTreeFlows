"""
TreeSBM benchmark suite (separate from core model code in `src/`).

Dependency direction is one-way: `benchmarks/` imports from `src/` (the core
model, TreeState, bridge modules), but `src/` never imports from `benchmarks/`.
Nothing here is on the training/inference path.

Layout:
  metrics/    model-independent metric library (topology, distributions,
              sequences, matching) — pure Python, unit-testable in isolation
  sim/        Track A synthetic simulators (custom forward BD + mutation)
  forecast/   Track B real-viral forecast-case construction
  tests/      benchmark unit tests
  (top level) generation harness + runnable track scripts (model-dependent)
"""
