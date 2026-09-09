# Group 003 — Gamma figure quantitative stats

Brazil geo held-out test group 3 (dates 2021-03-25 … 2021-05-10). Generated: `group_003_generated_matched.fasta` / `.nwk` (250 leaves). Observed: 300 leaves in `group_003_observed_anc_aa.fasta`.

EVEscape: official `evescape` column (`data/covid/spike_rbd_evescape_antibody_properties.csv` / `evescape_spike_rbd.pt`). **Higher (less negative) = more escape.** Percentile = % of all RBD substitutions with score ≤ this mut. Strain score = **sum** of per-mut EVEscape over RBD muts vs Wuhan ([evescape.org](https://evescape.org/) convention). Also report **mean** over scored RBD muts and over Gamma RBD-triad muts present.

## 1. Gamma reference

| Item | Value |
|------|-------|
| Definition | Consensus of observed leaves with full P.1 RBD triad **K417T + E484K + N501Y** |
| n triad leaves | **214** / 300 |
| Consensus muts vs Wuhan | L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y, D614G, H655Y, T1027I, V1176F |
| Designated tip (= consensus) | **MZ397166.1** (Hamming 0 to consensus; all 12 Gamma signature sites) |
| Obs prevalence | K417T 289/300; E484K 217/300; N501Y 292/300; D614G 296/300 |

Gamma signature: L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y, D614G, H655Y, T1027I, V1176F.

## 2. Min-Hamming generated leaf vs Gamma

Full-spike AA Hamming on length-matched sequences (L=1273).

| Metric | Value |
|--------|-------|
| **Min Hamming → Gamma consensus** | **289** |
| **Identity** | **0.7730** (77.30%) |
| Min Hamming → designated tip MZ397166.1 | 289 |
| Min Hamming → any triad obs tip | 288 (tip **MZ477792.1**) |
| Closest gen tip id | `root_child_0_child_0_child_1_child_1_child_1_child_1_child_0_child_1|leaf` |
| Gen Hamming distribution (→ consensus) | min=289, p25=313, **median=322**, mean=320.7, p75=330, max=351 |
| Tree-level mean best-id (obs→best gen) | ≈0.765 (matches `CLOSEST_GROUPS.md` 0.7661) |

**RBD-only** (sites 331–531): min Hamming **39** / 201 (identity 0.8060) on `root_child_1_child_1_child_0_child_0_child_1_child_1_child_1_child_1|leaf` — carries **K417T+N501Y**.

> Note: Generated spikes carry many extra substitutions outside the VOC signature, so absolute full-spike Hamming is large (~290). “Closest” = **best among 250 gen leaves**, not near-identity to P.1.

## 3. EVEscape ranking of the closest gen leaf

| Score | Closest (min full Hamming) | Rank / 250 | Best gen leaf |
|-------|----------------------------|------------|---------------|
| evescape **sum** (all scored RBD muts) | -107.329 | **#31** | -90.163 (Hamming 313) |
| evescape **mean** (all scored RBD muts) | -2.333 | #99 | — |
| evescape **sum of Gamma RBD triad muts present** | -1.435 (N501Y only) | **#1** (tied w/ other N501Y-only) | K417T+N501Y leaf: -2.542 |
| evescape mean of ≥90th-pct muts | -1.065 | #136 | — |

**Careful wording for claim 2:** Among gen leaves in the closest decile (Hamming ≤ 304), the **K417T+N501Y** leaf has the **highest evescape_sum** (overall rank **#5**, RBD Hamming **39** — best RBD match to Gamma). The absolute min-Hamming leaf is **not** #1 by full RBD evescape_sum (rank #31).

Gamma consensus (3 RBD muts): evescape_sum = -3.372, mean = -1.124.

## 4. High-escape / PMC / Gamma mutations recovered

### On min-Hamming gen leaf

Gamma / PMC sites present: **D138Y, N501Y, D614G, H655Y**

| Mutation | Present? | EVEscape | Percentile |
|----------|----------|----------|------------|
| E484K | no | -0.829 | 99.0 |
| K417T | no | -1.107 | 93.7 |
| N501Y | **yes** | -1.435 | 82.6 |
| D614G | **yes** | — (outside RBD) | — |
| D138Y | **yes** | — (outside RBD) | — |
| H655Y | **yes** | — (outside RBD) | — |

Additional ≥90th-pct EVEscape muts on this leaf: R346T (99.6%), Q493G (93.8%), K458N (91.8%), T470I (90.4%)

### On best-RBD / K417T+N501Y gen leaf

Tip: `root_child_1_child_1_child_0_child_0_child_1_child_1_child_1_child_1|leaf`
Gamma / PMC: **D138Y, K417T, N501Y, H655Y**
Full Hamming 303; RBD Hamming 39; evescape_sum rank #5.

≥90th-pct muts: K417T (93.7%), I472Q (90.7%)

### Gen-leaf recovery rates (n=250)

| Mut | n leaves | frac | EVEscape | pct | On closest? |
|-----|----------|------|----------|-----|-------------|
| L18F | 30 | 0.120 | — | — | N |
| T20N | 3 | 0.012 | — | — | N |
| P26S | 24 | 0.096 | — | — | N |
| D138Y | 177 | 0.708 | — | — | Y |
| R190S | 2 | 0.008 | — | — | N |
| K417T | 1 | 0.004 | -1.107 | 93.7 | N |
| E484K | 0 | 0.000 | -0.829 | 99.0 | N |
| N501Y | 220 | 0.880 | -1.435 | 82.6 | Y |
| D614G | 214 | 0.856 | — | — | Y |
| H655Y | 166 | 0.664 | — | — | Y |
| T1027I | 6 | 0.024 | — | — | N |
| V1176F | 94 | 0.376 | — | — | N |

**E484K was not recovered on any generated leaf (0/250).** K417T appears on **1** leaf (best RBD match). N501Y is nearly ubiquitous (220/250) and present at the gen **root** with D614G, D138Y, H655Y, V1176F.

### Trajectory (root → min-Hamming leaf): high-escape / Gamma / PMC acquisitions

| Mutation | First appears (node) | EVEscape | pct |
|----------|----------------------|----------|-----|
| D138Y | `root` | — | — |
| N501Y | `root` | -1.435 | 82.6 |
| D614G | `root` | — | — |
| H655Y | `root` | — | — |
| V1176F | `root` | — | — |
| K458R | `root_child_0_child_0` | -1.214 | 90.8 |
| R346T | `root_child_0_child_0_child_1` | -0.757 | 99.6 |
| Q498Y | `root_child_0_child_0_child_1` | -1.107 | 93.7 |
| Q493G | `root_child_0_child_0_child_1_child_1` | -1.105 | 93.8 |
| K458N | `root_child_0_child_0_child_1_child_1_child_1_child_1_child_0_child_1` | -1.173 | 91.8 |
| T470I | `root_child_0_child_0_child_1_child_1_child_1_child_1_child_0_child_1` | -1.226 | 90.4 |

## 5. Paste-ready caption numbers

- Gamma reference: consensus of **214** observed P.1-triad leaves; tip **MZ397166.1** (identical to consensus; full signature K417T/E484K/N501Y/D614G).
- Closest TreeSBM leaf: Hamming **289** / 1273 (**77.30%** identity) vs Gamma consensus; **288** vs nearest observed Gamma tip (**MZ477792.1**). Median gen Hamming **322**.
- Closest leaf recovers **N501Y** (EVEscape pct **82.6**), **D614G**, **D138Y**, **H655Y**; does **not** recover E484K or K417T.
- Best RBD match recovers **K417T+N501Y** (RBD Hamming **39**; evescape_sum rank **#5/250**).
- High-escape muts on closest-leaf path: **N501Y** (82.6%), **R346T** (99.6%), **Q493G** (93.8%), **K458N** (91.8%), **T470I** (90.4%).

## Files

- Per-leaf table: `group_003_gamma_gen_leaf_stats.csv`
- This summary: `group_003_gamma_figure_stats.md`
- Compact JSON: `group_003_gamma_figure_stats.json`
- **Trajectory comparison (best-leaf choice, obs↔gen mut intersection, caption):** `group_003_trajectory_compare.md`
