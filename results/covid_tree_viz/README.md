# COVID tree viz artifacts

Paired **observed** vs **model-generated** Newick trees for held-out COVID **group 40** (same root).
**Pair by leaf count** for side-by-side viz — do not compare full observed (207) against the fast 48-leaf gen.

| Field | Value |
|-------|-------|
| Group | 40 (`group_040`) |
| Split | `data/covid/test` (held-out test) |
| Checkpoint | `checkpoints/covid_v7_pmc_hotspot/best.pt` (epoch 36) |
| Mutation rate scale (mrs) | 0.5 |
| Fast gen | n_steps=40, max_leaves=48, seeds 42–44 (job 7437188) |
| Matched gen | n_steps=160, max_leaves=207, seed 42, `GEN_TAG=matched` (job 7437255) |
| branch_rate_scale | 6.0 |
| max_seq_len | 1280 |

## Recommended viz pairs

1. **Preferred (full size, both 207 tips):**  
   `group_040_observed.nwk` ↔ `group_040_generated_matched.nwk`
2. **Fast size-matched (both 48 tips):**  
   `group_040_observed_subsample.nwk` ↔ `group_040_generated.nwk`  
   Subsample = random tip prune of observed (seed=42), same root/group.

## Leaf counts

| File | Leaves |
|------|--------|
| `group_040_observed.nwk` | 207 |
| `group_040_generated_matched.nwk` | 207 |
| `group_040_generated_matched_s42.nwk` | 207 |
| `group_040_observed_subsample.nwk` | 48 |
| `group_040_generated.nwk` | 48 |
| `group_040_generated_s42.nwk` | 48 |
| `group_040_generated_s43.nwk` | 48 |
| `group_040_generated_s44.nwk` | 48 |

## Files

- `group_040_observed.nwk` — full ground-truth rooted tree
- `group_040_observed_subsample.nwk` — observed pruned to 48 random tips (seed 42)
- `group_040_generated.nwk` — fast primary sample (seed 42)
- `group_040_generated_s{SEED}.nwk` — additional fast samples
- `group_040_generated_matched.nwk` — leaf-count-matched generation (207 tips)
- optional `.fasta` sidecars with node AA sequences

## Quick viz

```bash
# Preferred full-size pair
python - <<'PY'
from ete3 import Tree
obs = Tree("group_040_observed.nwk", format=1)
gen = Tree("group_040_generated_matched.nwk", format=1)
print("leaves", len(obs), len(gen))
PY

# Fast size-matched pair
python - <<'PY'
from ete3 import Tree
obs = Tree("group_040_observed_subsample.nwk", format=1)
gen = Tree("group_040_generated.nwk", format=1)
print("leaves", len(obs), len(gen))
PY

# Or upload both .nwk files to https://itol.embl.de/
```

## Regenerate matched tree on Betty

```bash
GROUP=40 N_SAMPLES=1 MRS=0.5 MAX_LEAVES=207 N_STEPS=160 GEN_TAG=matched \
  SUBSAMPLE_N=48 REWRITE_README=1 \
  sbatch scripts/slurm_covid_tree_viz.sh
```

## Closest additional groups

Selected 2026-08-07 from held-out test (`covid_v7_pmc_hotspot`, mrs=0.5). See `CLOSEST_GROUPS.md`.

Prefer **equal-leaf** pairs (gen capped at 250):

- group 9 (closest): `group_009_observed_matched.nwk` ↔ `group_009_generated_matched.nwk` (combo=0.5304)
- group 36: `group_036_observed_matched.nwk` ↔ `group_036_generated_matched.nwk` (combo=0.5419)
- group 3: `group_003_observed_matched.nwk` ↔ `group_003_generated_matched.nwk` (combo=0.5627)

Group 40 kept as before (combo=0.8855 — much farther on terminal-edit / identity).
