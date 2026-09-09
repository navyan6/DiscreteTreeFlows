# Ancestral protein translation: which repair each virus needed

Three separate defects were corrupting `group_*_anc_aa.fasta`, and they need
three different fixes. Applying the SARS-CoV-2 repair everywhere would have made
flu and HIV worse, so the mode for each dataset is chosen from measurements, not
from the name.

Evidence: `scripts/panviral/audit_frame_offset.py` (is the alignment codon
framed?) and `scripts/panviral/probe_cds_start.py` (does the window open on the
start codon?). Assignments live in `MODE_BY_DATASET` in
`scripts/retranslate_anc_aa.py`.

## The measurements

Internal stop codons per sequence when translating in alignment coordinates. A
correct frame gives ~0; a wrong one scatters a stop every ~21 codons.

| dataset | width | %3 | offset 0 | offset 1 | offset 2 | verdict |
|---|---|---|---|---|---|---|
| covid (all splits) | 3822 | 0 | **0.00** | 54.98 | 144.47 | reference anchored |
| bdbv_temporal | 2700 | 0 | **0.00** | 56.00 | 67.00 | reference anchored |
| h3n2_epidemic | 1811 | 2 | 37.87 | 34.73 | 13.93 | no valid frame |
| h1n1_temporal | 1760 | 2 | 40.77 | 26.43 | 10.20 | no valid frame |
| hiv_geo | 2802 | 0 | 17.33 | 70.33 | 68.18 | no valid frame |

Where the CDS actually starts, against what `seq.find('ATG')` assumes:

| dataset | starts with ATG | 1st ATG | %3 | len@ATG | len@col0 |
|---|---|---|---|---|---|
| covid | 44/44 | 0 | 0 | 1273 | 1273 |
| bdbv | 0/42 | 113 | **2** | **24** | 905 |
| h3n2 | 0/44 | 51 | 0 | 566 | 603 |
| h1n1 | 40/44 | 0 | 0 | 566 | 578 |
| hiv | 44/44 | 0 | 0 | 871 | 881 |

## The three defects

**1. Frameshift from stripping gaps — SARS-CoV-2 only.**
The spike alignment is reference anchored, so three columns are one reference
codon. Stripping a gap run that is not a multiple of three shifts everything
downstream. Rare (0.1-0.7% of sequences) but concentrated in roots. Fixed by
translating column-wise: `mode="aligned"`.

**2. Truncation at an internal stop — flu and HIV.**
These alignments are per-group MAFFT nucleotide alignments. No column offset is
in frame for every sequence at once, and two widths are not even multiples of
three, so column-wise translation is meaningless there — it would manufacture
runs of `X`. Each sequence's own bases carry the frame, and stripping gaps
recovers the canonical protein: flu comes out at 566 aa, the textbook HA0
length. Their only bug was truncating at the first stop, which had cut HIV roots
to a fifth of their leaves' length. Fixed by `mode="ungapped"`, which keeps the
protein full length and marks internal stops `X`.

**3. Reading frame lost to the ATG search — Ebola.**
`nt_to_aa` locates the CDS with `seq.find('ATG')`, which is only right when the
window opens on the start codon. Bundibugyo's does not: its first ATG is at
index 110-113 and `113 % 3 == 2`, so every sequence was read out of frame and
died at the first in-frame stop. The result was 17-24 residue proteins from a
900-codon alignment containing no gaps at all — the gap logic was never even
involved. H3N2's first ATG is at 51, divisible by three, so it got the right
frame by luck. Fixed with `cds_start=0`.

## Result

49 dataset directories, 3694 trees. 717 truncated roots repaired. Re-running the
retranslation now reports zero changes everywhere, so the data is a fixed point.

`data/covid/test` was already repaired and reports 0 of 21001 sequences changed,
which is the check that the flu and HIV work did not disturb it.

Residual: 13 roots (0.35%) remain at 70-88% of leaf length, all in flu. Their
nucleotides carry real deletions, so this is an ancestral-reconstruction
artifact, not a translation one, and no translation mode recovers them.

Quarantined to `quarantine_malformed/` in `h1n1_epidemic/train`, all genuinely
broken rather than mistranslated:

- `group_052` — 468 nt for a gene whose modal width is 1701 (28%)
- `group_109` — root carries a 285 nt deletion; no mode reaches leaf length
- `group_016` — root translates to `MTMTHSINLLKYNHNRK...`, which is not HA at
  all (H1 HA begins `MKAILVVLLY`)

## Going forward

`stage_translate` in `scripts/run_all_groups.py` reads `TREESBM_TRANSLATE_MODE`
and `TREESBM_CDS_START`, and warns when the median protein falls below 80% of
the codons its alignment implies. That is the check that would have caught Ebola
immediately. New panviral viruses must be run through
`scripts/panviral/audit_frame_offset.py` and `probe_cds_start.py` before their
mode is set — `resolve_mode` refuses unknown datasets rather than guessing.
