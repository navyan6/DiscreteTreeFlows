# Closest COVID gen↔obs groups

Screened 10 held-out test groups (excl. 40) with `covid_v7_pmc_hotspot` @ mrs=0.5,
then size-matched regenerate for top-3 (max_leaves=min(obs,250), n_steps=160, 2 seeds).

**Combo distance** (lower = closer):
`0.45·RF + 0.35·terminal_edit + 0.20·(1−mean_best_identity)`
(+ 15% quartet when available). Observed tip-subsampled to gen leaf count
(seed=42) before sequence-matched RF / quartet / terminal-edit.

## Matched ranking (promoted)

| Rank | Group | Leaves (gen/obs) | RF ↓ | Quartet ↓ | Term-edit ↓ | Mean best-id ↑ | Cov@2% ↑ | Combo ↓ |
|-----:|------:|-----------------:|-----:|----------:|------------:|---------------:|---------:|--------:|
| 1 | 9 | 250/266 | 0.9490 | 0.6605 | 0.1753 | 0.9045 | 0.980 | 0.5304 |
| 2 | 36 | 250/291 | 0.8904 | 0.6300 | 0.2429 | 0.7964 | 0.000 | 0.5419 |
| 3 | 3 | 250/300 | 0.9505 | 0.5508 | 0.2580 | 0.7661 | 0.000 | 0.5627 |

## Recommended viz pairs (equal leaf count)

Because gen is capped at 250 tips, pair the **size-matched observed** Newick with the generated tree:

| Group | Observed (matched tips) | Generated |
|------:|-------------------------|-----------|
| 9 | `group_009_observed_matched.nwk` (250) | `group_009_generated_matched.nwk` (250) |
| 36 | `group_036_observed_matched.nwk` (250) | `group_036_generated_matched.nwk` (250) |
| 3 | `group_003_observed_matched.nwk` (250) | `group_003_generated_matched.nwk` (250) |

Full observed trees are also kept as `group_XXX_observed.nwk` (266 / 291 / 300 tips).

## Reference: group 40 (kept; not a closeness winner)

| Group | Leaves | RF | Quartet | Term-edit | Mean best-id | Cov@2% | Combo |
|------:|-------:|---:|--------:|----------:|-------------:|-------:|------:|
| 40 | 207/207 | 0.9896 | 0.6867 | 0.8694 | 0.1450 | 0.000 | 0.8855 |

Pair: `group_040_observed.nwk` ↔ `group_040_generated_matched.nwk`.

## Screen shortlist (96-tip, for reference)

See `SCREEN_RANKING.md`. Top screen groups were also 9 → 36 → 3 (same order as matched).

## Repro

Betty job `7437481`. Scripts: `scripts/slurm_covid_closest_groups.sh`,
`scripts/rank_covid_closest_groups.py`.
