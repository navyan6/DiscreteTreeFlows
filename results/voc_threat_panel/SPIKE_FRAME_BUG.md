# Spike deletions break reference-coordinate mutation calls

Found while diagnosing why `south_africa_beta` had 120 leaves labelled B.1.351
but only 3 apparently carrying K417N/E484K/N501Y.

## Cause

`scripts/covid_extract_spike.py` slices the nextclade-aligned genome at Spike's
reference coordinates and then strips gaps:

```13:13:scripts/covid_extract_spike.py
        region = str(rec.seq)[SPIKE_START:SPIKE_END].replace("-", "").upper()
```

Stripping the gaps is what breaks the frame. A lineage with a Spike deletion
becomes shorter, so every residue downstream of the deletion sits at a lower
index than its reference position. Anything that calls mutations by reference
position — `voc_threat_lib.muts_present`, the screens, the leakage audit —
silently misses them.

## Evidence

`scripts/check_spike_frame_shift.py`, counting target-lineage leaves carrying
each defining mutation at the reference index vs shifted indices:

| case | mutation | at ref | at ref−3 |
|---|---|---|---|
| Beta / B.1.351 | K417N | 3 | **114** |
| Beta / B.1.351 | E484K | 1 | **109** |
| Beta / B.1.351 | N501Y | 3 | **110** |
| Alpha / B.1.1.7 | N501Y | 2 | **118** |
| Alpha / B.1.1.7 | P681H | 2 | **118** |
| Epsilon / B.1.427/429 | L452R | **118** | 0 |
| Epsilon / B.1.427/429 | W152C | **117** | 0 |

Alpha deletes ΔH69/V70 and ΔY144 (3 residues) and Beta deletes Δ242–244 (3
residues), both upstream of their defining sites — hence the uniform −3. Epsilon
has no Spike deletion and lines up exactly, which is the control.

## Blast radius

Affects any lineage with a Spike deletion, everywhere reference positions are
used:

| VOC | Spike deletion | affected |
|---|---|---|
| Alpha | Δ69–70, Δ144 | yes (−3) |
| Beta | Δ242–244 | yes (−3) |
| Delta | Δ156–157 | yes (−2, so L452R and T478K both) |
| Omicron BA.1 | Δ69–70, Δ143–145, Δ211, ins214EPE | yes, net shift varies |
| **Gamma** | none | **no** |
| Epsilon, Iota, Mu | none / insertion only | no |

So the Gamma Brazil case study — the current headline — is unaffected. But these
are undercounts and need recomputing:

- `train_leakage_audit.json` counts for Delta (1321), Omicron (774), Alpha (870)
- `bundle_frac` for the Brazil Delta cases g007 / g011
- the Alpha and Beta rows of the new-roots panel, which looked unusable and are
  in fact fine

## Fix

Correct this in the evaluation layer, not the extraction. The model can keep
operating on ungapped sequences; only our mutation calling needs reference
coordinates. Pairwise-align each sequence to the 1273-aa Spike reference and map
reference positions through the alignment, cached by sequence hash — there are
only ~10,001 distinct Spike strings across all 146,495 training sequences, so
this is cheap.

Re-extracting with gaps preserved would instead change the model's inputs and
force a retrain, for no modelling benefit.

Not a bug: `data/covid/wt.txt` is a FASTA file, so its raw character count is
1410 while the Spike protein it holds is the correct 1273 aa.
