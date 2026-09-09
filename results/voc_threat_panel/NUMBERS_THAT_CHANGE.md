# Every number that moves after the Spike frame fix

Cause and evidence: `SPIKE_FRAME_BUG.md`. Fix: `voc_threat_lib.ref_pos_map` /
`seq_at` now map reference positions through a pairwise alignment to the 1273-aa
Spike reference, so deletion-carrying lineages are read at the right residue.
Sequences of length 1273 skip alignment entirely, so the fix is nearly free.

**Gamma is the control throughout.** P.1 has no Spike deletion, so if the fix were
wrong Gamma would move. It does not.

---

## 1. Confirmed changed — recomputed

### Training-set VOC carriers (`train_leakage_audit.json` → `_fixed.json`)

335 trees, 81,271 leaves.


| VOC          | trees before → after | carrier leaves before → after | factor          | frac of train leaves |
| ------------ | -------------------- | ----------------------------- | --------------- | -------------------- |
| **Gamma**    | 36 → 37              | 232 → **241**                 | 1.04× (control) | 0.003                |
| Beta         | 4 → 85               | 54 → **3,568**                | 66×             | 0.044                |
| Delta        | 46 → 244             | 1,321 → **38,550**            | 29×             | **0.474**            |
| Omicron BA.1 | 20 → 180             | 774 → **30,756**              | 40×             | **0.378**            |
| Alpha        | 38 → 212             | 870 → **29,626**              | 34×             | **0.365**            |
| Lambda       | 3 → 6                | 3 → **55**                    | 18×             | 0.001                |
| Mu           | 48 → 179             | 330 → **10,240**              | 31×             | 0.126                |
| Zeta         | 69 → 142             | 621 → **4,564**               | 7×              | 0.056                |


The framing consequence is larger than the arithmetic: Delta's signature is
present in **47% of all training leaves**, Omicron's in 38%, Alpha's in 36%. Any
sentence implying these signatures were rare in training is wrong. Under the
tree-level standard this is still acceptable — shared mutations are expected —
but it must be stated accurately.

### Brazil held-out cases (paper-facing)


| case            | metric                          | before → after                           |
| --------------- | ------------------------------- | ---------------------------------------- |
| **Gamma g003**  | `bundle_frac`                   | 0.7133 → **0.7167** (control, unchanged) |
| Gamma g003      | `exact_mut_recall`              | 1.000 → 1.000                            |
| Gamma g003      | root bundle / prospective       | `['N501Y']` / True — unchanged           |
| Delta g007      | `bundle_frac`                   | 0.0100 → **0.0533** (5.3×)               |
| Delta g011      | `bundle_frac`                   | 0.0334 → **0.4181** (12.5×)              |
| Delta g007/g011 | `exact_mut_recall`, prospective | 1.000, True — unchanged                  |


### New unseen-root panel (`voc_roots_seq_overlap.json`)

`check_voc_roots_overlap.py` was still doing raw `s[pos-1]` indexing and so
bypassed the fix; it now calls `voc_threat_lib.seq_at`. Regenerated values:


| population        | tree      | bundle-carrying leaves before → after            |
| ----------------- | --------- | ------------------------------------------------ |
| south_africa_beta | g001      | 3 → **143**                                      |
| uk_alpha          | g001      | 2 → **123**                                      |
| uk_alpha          | g002      | 3 → **124**                                      |
| uk_alpha          | g003      | 6 → **123**                                      |
| india_delta       | g001      | 0 → 1                                            |
| usa_epsilon       | g001–g003 | 116/117/117 → 118/118/120 (no deletion, control) |
| usa_iota          | g001–g003 | 68/78/79 → 69/79/80 (no deletion, control)       |
| colombia_mu       | g001      | 124 → 125                                        |


Alpha and Beta looked unusable and are in fact among the best-populated trees.
Epsilon and Iota barely move, which is the expected control behaviour — neither
lineage carries a Spike deletion.

Note the second-order consequence: most of the newly-found Alpha and Beta bundle
leaves are byte-identical to a training Spike (only 15–26 novel per tree), so
these trees gain population but not prospective strength. Iota and Mu remain the
Tier-1 cases at 100% novel. See `README.md` for the tiering.

---

## 2. Confirmed *not* changed

- **The Gamma Brazil case study**, i.e. the current headline figure and
`group_003_gamma_figure_stats.md` / `group_003_trajectory_compare.md`. P.1 has
no Spike deletion; recomputed values match to within one leaf.
- Any metric that never indexes by reference position: Hamming/edit-distance
coverage curves, tree-topology metrics, PLL scores. These compare sequences to
each other, not to reference coordinates.
- `data/covid/wt.txt` — not a bug; it is a FASTA file, so its raw character
count (1410) differs from the 1273-aa protein it holds.

---

## 3. Answered: is column-wise leaf-vs-root comparison safe?

**No, not in aggregate — but the trees carrying the paper's claims are fine.**

`check_tree_frame_consistency.py` compares every leaf to its tree's root. Equal
length is necessary but not sufficient (two sequences can both be 1271 aa via
different deletions), so it also requires Hamming ≤ 50 against the root; a true
frame shift smears the downstream sequence into hundreds of mismatches.


| split         | leaves | share the root's frame |
| ------------- | ------ | ---------------------- |
| train         | 81,271 | 38,762 (**47.7%**)     |
| test (Brazil) | 11,701 | 5,463 (**46.7%**)      |


The per-tree distribution is bimodal rather than middling — trees are almost
always all-consistent or all-broken:


| frac of leaves sharing root frame | train trees |
| --------------------------------- | ----------- |
| 0.0                               | **130**     |
| 0.1–0.9                           | 95          |
| 1.0                               | **110**     |


The 130 zero-consistency trees are the expected biology: the inferred ancestral
root carries no deletion while every descendant does, so root and leaves sit in
different frames by construction. Column-wise comparison in those trees is
meaningless.

Among equal-length pairs the frame is almost always genuinely shared — median
Hamming to root is 2 (train) / 3 (test), and only 730 of 39,492 equal-length
pairs exceed the shift threshold.

### The three paper-facing Brazil trees


| tree                      | leaves | root len | frac sharing root frame | median Hamming |
| ------------------------- | ------ | -------- | ----------------------- | -------------- |
| **Gamma g003 (Figure 3)** | 300    | 1273     | **0.987**               | 7              |
| Delta g007                | 300    | 1273     | 0.953                   | 8              |
| Delta g011                | 299    | 1273     | 0.595                   | 9              |


So Figure 3 is safe on this axis too, and Delta g007 nearly so. **Delta g011 is
the one to treat with caution** — 40% of its leaves are read in the wrong frame
by any column-wise metric.

### What this means

- Any *aggregate* COVID metric computed leaf-vs-root by column — the site
recall, `aa_acc | hit`, and EVEscape enrichment in
`eval_evescape_enrichment.py` — is computed on misaligned columns for roughly
half its inputs, and for all inputs in 39% of trees.
- Per-case figures built on Gamma g003 or Delta g007 are unaffected in practice.

## 4. Still open, NOT yet quantified

**Reference-coordinate masks.** The RBD hotspot mask, the curated PMC literature
mask, and the EVEscape score tensor are all indexed by Spike position and applied
to padded columns. Within-tree frame consistency does not help here: the mask is
defined in reference coordinates, so for the **74% of sequences that are not
1273 aa** it lands on the wrong residues regardless of how internally consistent
the tree is.


| split         | leaves | length 1273  | length != 1273   |
| ------------- | ------ | ------------ | ---------------- |
| train         | 81,271 | 21,282 (26%) | **59,989 (74%)** |
| test (Brazil) | 11,701 | 3,859 (33%)  | **7,842 (67%)**  |


Most common train lengths: 1271 (26,454), 1273 (21,282), 1268 (10,767),
1270 (9,733), 1267 (8,617).

This would touch the EVEscape-guided ablations rather than only the reporting,
since the masks feed training as well as evaluation. Quantifying it means
re-scoring with `ref_pos_map` applied to the mask lookup, which has not been run.

---

## 5. How far this reaches beyond COVID

`scripts/audit_frame_all_datasets.py` → `frame_audit_all_datasets.json`.

Because TreeSBM only substitutes, a generated leaf always inherits the root's
frame. So the two numbers that matter per dataset are how often an *observed*
leaf shares its root's frame (decides column-wise leaf metrics) and how often the
*root* sits at reference length (decides reference-coordinate masks).


| Dataset                        | trees | modal len | roots at modal | leaves at modal | leaves in root's frame |
| ------------------------------ | ----- | --------- | -------------- | --------------- | ---------------------- |
| SARS-CoV-2 Spike (Brazil test) | 40    | 1271      | 0.225          | 0.465           | **0.467**              |
| SARS-CoV-2 Spike (train)       | 335   | 1271      | 0.164          | 0.326           | **0.477**              |
| Influenza H3N2 HA              | 45    | 566       | 0.889          | 0.914           | **0.936**              |
| Influenza H1N1 HA              | 40    | 566       | 0.725          | 0.872           | **0.744**              |
| HIV Env (geo test)             | 14    | 853       | 0.071          | 0.087           | **0.018**              |


**HIV Env is the worst case by a wide margin.** Env's variable loops differ in
length between essentially every isolate, so only 1.8% of leaves share their
root's frame — every column-wise HIV number currently in the paper is noise.
H3N2 is nearly immune at 94%. H1N1 is real but moderate.

The paper-level consequences are worked through in
`[PAPER_NUMBERS_TO_UPDATE.md](PAPER_NUMBERS_TO_UPDATE.md)`, which lists the
affected cells table by table and the CPU-only rescoring path that fixes them.

---

## 6. The fix, verified

`benchmarks/metrics/sequences.py::positional_recovery` now maps the observed and
generated leaves onto the root's frame by pairwise alignment before comparing;
equal-length sequences take the identity map and skip alignment entirely.
`any_descendant_mut_recovery` got the same treatment. `align=False` reproduces
the old behaviour for back-comparison.

`scripts/test_positional_recovery_frame.py` — all checks pass, 59/59 existing
benchmark tests still pass.

**Synthetic control** (deletion upstream of a known substitution, generated leaf
correct at the homologous residue): `mut_recovery` 0.000 → **1.000**. The legacy
metric scored a perfect forecast as a total miss.

**Gamma g003** — 147 of 150 pairs are in-frame, and on those the two code paths
are **bit-identical** (0 mismatches). Aggregates barely move:


| metric             | legacy → aligned |
| ------------------ | ---------------- |
| `mut_recovery`     | 0.0147 → 0.0150  |
| `cons_retention`   | 0.7788 → 0.7768  |
| `aa_acc_given_hit` | 0.0768 → 0.0754  |


**Off-frame contrast, g040** (2% of leaves in frame) — large movement, and it is
**downward**:


| metric             | legacy → aligned    |
| ------------------ | ------------------- |
| `mut_recovery`     | 0.0760 → **0.0157** |
| `cons_retention`   | 0.8468 → **0.7111** |
| `aa_acc_given_hit` | 0.2524 → **0.0411** |


This is the direction to expect for the paper. Legacy scoring inflates these
metrics on off-frame trees, because a frame shift makes almost every position
look like a "mutating site" and the model then collects chance-level credit
across thousands of spurious sites. Removing them shrinks the denominator to the
handful of real mutations, where the model does much worse. **Table 4's
SARS-CoV-2 `AA|hit` (0.375) and key mutation recall (0.219) should be expected to
fall, not rise.**

Two trees is an indication of direction and rough magnitude, not a forecast of
the final numbers.

---

## 7. Found while testing: truncated ancestral roots

28% of COVID trees and 64% of HIV trees have a root that is a short N-terminal
fragment (as little as 6–27 aa) while their leaves are full length. Generated
leaves inherit the root's length, so those trees can never contain the RBD at
all. This is a data defect, it is larger than the frame bug, and it also reaches
the distance-based metrics the frame fix leaves alone.

See `[TRUNCATED_ROOTS.md](TRUNCATED_ROOTS.md)`. It should be resolved *before*
the rescoring pass.