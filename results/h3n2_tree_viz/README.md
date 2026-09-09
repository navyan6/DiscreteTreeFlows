# H3N2 tree viz artifacts

Paired **observed_matched** vs **generated_matched** Newick trees for 2 held-out H3N2 test groups.
Same root/group; observed tip-subsampled to generated leaf count (seed=42) for side-by-side viz.

| Field | Value |
|-------|-------|
| Split | `data/h3n2/test` (held-out test) |
| Checkpoint | `checkpoints/h3n2_v3_lit_hotspot/best.pt` (epoch 57) |
| Mutation rate scale (mrs) | 0.5 |
| Match max_leaves / n_steps | 250 / 160 |
| Screen candidates | 46, 5, 20, 40, 55, 59 (96 tips / 60 steps) |
| Winners (closest screen combo) | **55**, **40** |
| branch_rate_scale | 6.0 |
| max_seq_len | 566 |
| Job | SLURM 7439927 |

## Selected groups

| Group | Dates | Full obs leaves | Matched leaves | Combo ↓ | Viz pair |
|------:|-------|----------------:|---------------:|--------:|----------|
| 55 | 2024-11-02..2024-11-09 | 400 | 250/250 | 0.5433 (screen) / 0.5452 (matched) | `group_055_observed_matched.nwk` ↔ `group_055_generated_matched.nwk` |
| 40 | 2024-06-21..2024-06-28 | 399 | 250/250 | 0.5451 (screen) / 0.5450 (matched) | `group_040_observed_matched.nwk` ↔ `group_040_generated_matched.nwk` |

Location: not in CSV (name/date only); GISAID H3N2 HA sequences.

## Recommended viz pairs (equal leaf count)

1. **group 55:** `group_055_observed_matched.nwk` ↔ `group_055_generated_matched.nwk` (both 250 tips)
2. **group 40:** `group_040_observed_matched.nwk` ↔ `group_040_generated_matched.nwk` (both 250 tips)

Full observed trees kept as `group_XXX_observed.nwk` (399 / 400 tips).

## Leaf counts

| File | Leaves |
|------|--------|
| `group_040_observed.nwk` | 399 |
| `group_040_observed_matched.nwk` | 250 |
| `group_040_generated_matched.nwk` | 250 |
| `group_055_observed.nwk` | 400 |
| `group_055_observed_matched.nwk` | 250 |
| `group_055_generated_matched.nwk` | 250 |

## Quick viz

```bash
python - <<'PY'
from ete3 import Tree
for g in ("055", "040"):
    obs = Tree(f"group_{g}_observed_matched.nwk", format=1)
    gen = Tree(f"group_{g}_generated_matched.nwk", format=1)
    print(g, "leaves", len(obs), len(gen))
PY
# Or upload both .nwk files to https://itol.embl.de/
```

## Regenerate on Betty

```bash
sbatch scripts/slurm_h3n2_tree_viz.sh
```
