# Eval metrics: where they live + leaf vs tree

## Primary KPIs (leaf-only)

| Metric | Definition | Scope |
|--------|------------|-------|
| `mut_recovery` | P(gen==GT \| root≠GT) | **Leaf**: best-match gen leaf vs GT leaf vs root |
| `site_recall` | P(gen≠root \| root≠GT) | same |
| `aa_acc_given_hit` | P(gen==GT \| root≠GT & gen≠root) | same |
| `cons_retention` | P(gen==root \| root==GT) | same |
| `pmc_hotspot_mut_frac` / `flu_hotspot_mut_frac` | \|mut_cols ∩ lit\| / \|mut_cols\| | **Leaf** gen≠root / GT≠root |
| `ntd/rbd/rbm/s1/s2_mut_frac` | domain-band mut frac (COVID) | **Leaf** (same) |

Internal / ancestral nodes are **not** used in these primaries. That is intentional for tip forecasting: only the end-state leaf AA matters. If GT mutations appear only on terminal edges while the model mutates earlier (or vice versa), leaf scoring remains the fair check of recovered tip genotype.

## Where computed / called

| Location | Role |
|----------|------|
| `benchmarks/metrics/sequences.py` → `positional_recovery` | Canonical definition |
| `scripts/eval_single_tree.py` | Same formula (local copy) + GT leaf comparison |
| `scripts/eval_evescape_enrichment.py` | Aggregate COVID/flu eval + lit/domain fracs |
| `scripts/eval_leaf_holdout.py` | Held-out leaf recovery |
| `scripts/eval_test_set.py` | Gen smoke / identity (not full siteaa) |
| `benchmarks/track_a.py` | Track A: leaf pool + `positional_recovery` |
| `scripts/slurm_eval_covid*.sh`, `slurm_inference_sweep.sh`, … | Call enrichment |

## Extra tree-wide companions (do not replace primary)

Emitted by `eval_evescape_enrichment.py` and available in `benchmarks/metrics/sequences.py`:

- `mut_recovery_any_descendant` / `site_recall_any_descendant`
- `mut_site_recall_path_union` (via `path_union_mutation_recovery` / `mutation_pr_f1`)

Keep reporting **both** leaf primary and these extras when discussing fairness.
