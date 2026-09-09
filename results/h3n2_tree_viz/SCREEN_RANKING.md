# H3N2 screen ranking (96-tip)

Ranked by combo distance (lower = closer):
`0.45·RF + 0.35·terminal_edit + 0.20·(1−mean_best_identity)`
(+ 15% quartet when available). Observed tip-subsampled to gen leaf count
(seed=42) before matching.

| Rank | Group | Leaves (gen/obs) | RF ↓ | Quartet ↓ | Term-edit ↓ | Mean best-id ↑ | Cov@2% ↑ | Combo ↓ | Tag |
|-----:|------:|-----------------:|-----:|----------:|------------:|---------------:|---------:|--------:|-----|
| 1 | 55 | 96/400 | 0.9876 | 0.7405 | 0.1247 | 0.8980 | 0.000 | 0.5433 | screen |
| 2 | 40 | 96/399 | 0.9869 | 0.7408 | 0.1286 | 0.8928 | 0.000 | 0.5451 | screen |
| 3 | 46 | 96/378 | 0.9874 | 0.7235 | 0.1395 | 0.8851 | 0.000 | 0.5473 | screen |
| 4 | 5 | 96/396 | 1.0000 | 0.7877 | 0.1225 | 0.9015 | 0.000 | 0.5538 | screen |
| 5 | 59 | 96/400 | 1.0000 | 0.7133 | 0.1485 | 0.8765 | 0.000 | 0.5547 | screen |
| 6 | 20 | 96/399 | 1.0000 | 0.7883 | 0.1446 | 0.8774 | 0.000 | 0.5646 | screen |

