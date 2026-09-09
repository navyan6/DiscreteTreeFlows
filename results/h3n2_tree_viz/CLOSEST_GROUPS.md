# Closest H3N2 gen↔obs groups

Ranked by combo distance (lower = closer):
`0.45·RF + 0.35·terminal_edit + 0.20·(1−mean_best_identity)`
(+ 15% quartet when available). Observed tip-subsampled to gen leaf count
(seed=42) before matching.

| Rank | Group | Leaves (gen/obs) | RF ↓ | Quartet ↓ | Term-edit ↓ | Mean best-id ↑ | Cov@2% ↑ | Combo ↓ | Tag |
|-----:|------:|-----------------:|-----:|----------:|------------:|---------------:|---------:|--------:|-----|
| 1 | 40 | 250/399 | 1.0000 | 0.7370 | 0.1240 | 0.9112 | 0.000 | 0.5450 | matched |
| 2 | 55 | 250/400 | 1.0000 | 0.7316 | 0.1236 | 0.9045 | 0.000 | 0.5452 | matched |


## Promoted for viz

- group 40: `group_040_observed_matched.nwk` ↔ `group_040_generated_matched.nwk` (combo=0.5450, leaves=250)
- group 55: `group_055_observed_matched.nwk` ↔ `group_055_generated_matched.nwk` (combo=0.5452, leaves=250)
