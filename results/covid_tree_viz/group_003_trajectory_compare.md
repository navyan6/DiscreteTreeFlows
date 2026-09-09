# Group 003 — Trajectory comparison (obs vs gen) for Gamma figure

Companion to `group_003_gamma_figure_stats.md`. Goal: pick the most figure-convincing generated leaf, build root→leaf mutation tables against the observed Gamma tip, score shared muts with EVEscape, and recommend the narrative.

**Wuhan / roots**

| Sequence | Source | Gamma-signature muts vs Wuhan |
|----------|--------|-------------------------------|
| Wuhan-Hu-1 | `data/covid/wt.txt` (P0DTC2, L=1273) | — |
| Obs root | `NODE_0000001` in `group_003_observed_anc_aa.fasta` | D138Y, N501Y, D614G, H655Y, V1176F |
| Gen root | `root\|root` in `group_003_generated_matched.fasta` | D138Y, N501Y, D614G, H655Y, V1176F |
| Obs Gamma tip | **MZ397166.1** (= Gamma consensus; 12/12 signature) | full P.1 set |
| Trees | obs path: `group_003_observed.nwk` + anc fasta; gen path: leaf-id path in matched fasta/nwk |

Both trees already start from the same Brazil-background signature block. Distinctive Gamma acquisitions on the observed path are mainly NTD + RBD triad completions (P26S → K417T → … → **E484K**).

EVEscape: official `evescape` column; higher (less negative) = more escape. Percentile = % of all RBD substitutions with score ≤ this mut.

---

## 1. Recommendation (TL;DR)

**Primary gen leaf for the figure:**  
`root_child_1_child_1_child_0_child_0_child_1_child_1_child_1_child_1|leaf`  
(= min RBD Hamming to Gamma / only K417T+N501Y leaf)

**Matched observed tip:** `MZ397166.1` (Gamma consensus tip).

**Shared mutset (leaf − Wuhan):** **D138Y, K417T, N501Y, H655Y**  
Scored shared escape: **K417T (93.7th pct), N501Y (82.6th pct)**.

**Is shared-mut + EVEscape the strongest story?**  
**Yes for a single tip-vs-tip figure panel** — but only with this leaf, and only if the caption is honest that (i) E484K is never generated (0/250), (ii) D614G/V1176F are at the gen root then *lost* on this leaf (D614Y / V1176I), and (iii) three of the four shared muts are already at both roots; the distinctive de novo recovery is **K417T**.

**Stronger population-level backup (use in text / supplement, not instead of the table):** among all 250 gen leaves, **11/12** exact Gamma-tip substitutions appear at least once (only **E484K** missing); **12/12** Gamma sites are hit by some substitution. Group-3 mutrec (`track_a_covid_v5_mutrec_K100.json`): mut_recall **0.887**, site_recall **1.0**.

---

## 2. Best-leaf criteria bake-off

| Criterion | Tip (short) | Full Hamm | RBD Hamm | n Γ-sig /12 | Γ-sig hits | evescape_sum (rank/250) | Shared w/ MZ397166 | Verdict |
|-----------|-------------|-----------|----------|-------------|------------|-------------------------|--------------------|---------|
| Min full Hamming | `…child_0_child_1\|leaf` | **289** | 47 | 4 | D138Y N501Y D614G H655Y | −107.3 (**#31**) | D138Y N501Y D614G H655Y | Closest overall, but **no K417T**; weakest antigenic story |
| **Min RBD Hamming** | `…child_1_child_1\|leaf` | 303 | **39** | 4 | D138Y **K417T** N501Y H655Y | −95.0 (**#5**) | D138Y **K417T** N501Y H655Y | **PRIMARY** — only leaf with K417T; best RBD match; high escape on shared RBD muts |
| Max Γ-sig overlap | `…child_0_child_0\|leaf` | 292 | 44 | **6** | L18F D138Y N501Y D614G H655Y V1176F | −95.5 (#7) | L18F D138Y N501Y D614G H655Y V1176F | Best signature *count*, but **no RBD triad beyond N501Y** |
| Max evescape_sum | `…child_1_child_1\|leaf` (other clade) | 313 | 41 | 5 | L18F N501Y D614G H655Y V1176F | **−90.2 (#1)** | L18F N501Y D614G H655Y V1176F | Highest escape sum via gen-only RBD muts (L455R etc.), not Gamma triad |
| Tradeoff score* | = min-RBD leaf | — | — | — | — | — | — | **Winner** (core×3 + sig×0.5 + K417T − RBD/50 − rank/100) |

\*Ad hoc score used only to break ties; K417T leaf wins clearly (≈8.2 vs ≈5.1 for max-sig).

**Why not min full Hamming?** It recovers early background (incl. D614G) but misses the defining RBD change K417T that carries high EVEscape (93.7th pct). Absolute Hamming is dominated by ~290 non-Gamma substitutions on every gen leaf — a bad figure metric.

**Why not max Γ-sig?** Six shared signature sites look great in a checklist, but without K417T/E484K the EVEscape trajectory story collapses to N501Y alone.

---

## 3. Trajectory tables (paste-ready)

### 3A. Gamma-signature path: Observed MZ397166.1

Ordered by first appearance on `NODE_0000001` → … → tip (from Newick + anc AA).

| Order | Mutation | First node | In gen primary leaf? | EVEscape | pct |
|------:|----------|------------|----------------------|----------|-----|
| 1 | D138Y | root (`NODE_0000001`) | **yes** | — (NTD) | — |
| 2 | N501Y | root | **yes** | −1.435 | **82.6** |
| 3 | D614G | root | no (leaf has D614Y; present at gen root) | — | — |
| 4 | H655Y | root | **yes** | — | — |
| 5 | V1176F | root | no (leaf has V1176I; present at gen root) | — | — |
| 6 | P26S | `NODE_0000087` | no | — | — |
| 7 | **K417T** | `NODE_0000087` | **yes** | −1.107 | **93.7** |
| 8 | T1027I | `NODE_0000087` | no | — | — |
| 9 | T20N | `NODE_0000085` | no | — | — |
| 10 | R190S | `NODE_0000074` | no | — | — |
| 11 | L18F | `NODE_0000073` | no | — | — |
| 12 | **E484K** | `NODE_0000064` | **never in any gen leaf** | −0.829 | **99.0** |

### 3B. Gamma / high-escape path: Generated primary leaf (retained to tip)

Path reconstructed from leaf-id (`root` → `root_child_1` → … → leaf). Rows = mutations **retained on the leaf** that are Gamma-signature or ≥90th-pct EVEscape.

| Order | Mutation | First appears | Class | EVEscape | pct |
|------:|----------|---------------|-------|----------|-----|
| 1 | D138Y | `root\|root` | Γ-sig (shared) | — | — |
| 2 | N501Y | `root\|root` | Γ-sig / RBD (shared) | −1.435 | **82.6** |
| 3 | D614G | `root\|root` | Γ-sig; **reverted** before tip → D614Y | — | — |
| 4 | H655Y | `root\|root` | Γ-sig (shared) | — | — |
| 5 | V1176F | `root\|root` | Γ-sig; **reverted** → V1176I | — | — |
| 6 | I472Q | `root_child_1_child_1\|internal` | gen-only ≥90th pct | −1.216 | **90.7** |
| 7 | **K417T** | leaf (terminal step) | Γ-sig / RBD (shared) | −1.107 | **93.7** |

(L455R appears mid-path then reverts; not retained.)

### 3C. Intersection summary (obs tip ∩ gen primary leaf)

| Mutation | Shared? | EVEscape | pct | Note |
|----------|---------|----------|-----|------|
| D138Y | **yes** | — | — | Already at both roots |
| K417T | **yes** | −1.107 | **93.7** | **Best figure mut** — acquired on gen leaf; early on obs path |
| N501Y | **yes** | −1.435 | **82.6** | Already at both roots |
| H655Y | **yes** | — | — | Already at both roots |
| E484K | **no** | −0.829 | **99.0** | Missing from all 250 gen leaves |
| D614G | **no*** | — | — | *At gen root, lost on this leaf |
| L18F T20N P26S R190S T1027I V1176F | no | — | — | On max-sig leaf: L18F+V1176F+D614G present |

**Gen-only high-escape (primary leaf):** I472Q (90.7th pct).  
**Obs-only (missing) headline mut:** E484K (99.0th pct).

### 3D. Optional secondary leaf (max signature overlap)

Tip: `root_child_1_child_0_child_1_child_1_child_1_child_0_child_1_child_0_child_0|leaf`  
Shared with MZ397166: **L18F, D138Y, N501Y, D614G, H655Y, V1176F** (6/12).  
Scored shared escape: **N501Y only**. Use only if the figure prioritizes signature checklist over antigenic RBD.

---

## 4. Narrative recommendation (honest)

| Framing | Strength | Weakness | Use? |
|---------|----------|----------|------|
| **Shared-mut + EVEscape trajectory table** (primary leaf) | Concrete tip-level story; K417T+N501Y are real Gamma antigenic muts with high escape percentiles; path-reconstructable | Only 4 shared muts; 3/4 are root-shared; no E484K; D614G not retained | **Yes — main figure panel** |
| Set-level signature recovery (any of 250 leaves) | 11/12 exact tip muts; 12/12 sites; clear “almost full P.1 repertoire” | Not a single trajectory; K417T is rare (1/250) | Strong caption / results sentence |
| Min full-Hamming leaf | Lowest Hamm (289) | Misses K417T; EVEscape rank only #31 | Avoid as hero tip |
| Max evescape_sum leaf | Rank #1 escape | Escape driven by non-Gamma muts (L455R, G446D, …) | Misleading for “recovers Gamma” |
| Mutrec enrichment (group 3) | mut_recall 0.887, site_recall 1.0 | Aggregate metric; less visual | Methods / supplement |
| Path-union site recall | All Gamma sites touched somewhere | AA identity at site often wrong; E484K AA never hit | Don’t lead with this alone |

**Bottom line:** Lead with the **K417T+N501Y leaf trajectory table + EVEscape on shared muts**. Do **not** claim full Gamma triad recovery. Pair with one sentence on set-level 11/12 signature recall so the missing E484K reads as a specific failure mode, not total antigenic miss.

---

## 5. Caption (paste-ready)

> In Brazil held-out group 3, TreeSBM’s best RBD match to the P.1 (Gamma) consensus tip (MZ397166.1; RBD Hamming 39/201; evescape_sum rank 5/250) recovers the high-escape RBD substitutions **K417T** (EVEscape 93.7th percentile) and **N501Y** (82.6th), together with shared early Spike changes D138Y and H655Y. Side-by-side root→leaf tables show these mutations on both the observed Gamma path and the generated path (K417T acquired at the generated tip). The third triad mutation, **E484K** (99.0th percentile), is absent from all 250 generated leaves; D614G is present at the generated root but not retained on this tip. Across the full generated set, 11/12 exact Gamma-tip substitutions appear at least once (site recall 12/12).

---

## 6. Files

| File | Role |
|------|------|
| `group_003_generated_matched.fasta` / `.nwk` | Gen tree + ancestral AA |
| `group_003_observed_anc_aa.fasta` + `group_003_observed.nwk` | Obs path reconstruction |
| `group_003_gamma_gen_leaf_stats.csv` | Per-leaf Hamming / EVEscape / signature hits |
| `group_003_gamma_figure_stats.md` | Prior quantitative summary |
| This file | Leaf choice + trajectory intersection for the figure |
