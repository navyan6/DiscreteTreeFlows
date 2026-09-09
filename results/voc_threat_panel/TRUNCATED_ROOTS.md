# Truncated ancestral roots

Found while regression-testing the alignment fix. This is a **separate and larger
defect than the frame shift**, and it is not a metric bug — it is in the data.

## What it is

In a large fraction of trees, the reconstructed ancestral root is a short
N-terminal fragment of the protein rather than the full sequence, while the
leaves of the same tree are full length.

| Split | trees | roots >10% off the leaf length | leaves in those trees | shortest roots seen |
|---|---|---|---|---|
| SARS-CoV-2 Spike (train) | 335 | **94 (28.1%)** | 22,631 | 27, 28, 35, 46, 65 aa |
| SARS-CoV-2 Spike (test) | 40 | **11 (27.5%)** | 3,246 | 75, 98, 151, 152, 156 aa |
| HIV Env (geo test) | 14 | **9 (64.3%)** | 2,297 | 6, 26, 83, 131, 145 aa |
| Influenza H1N1 HA | 40 | **8 (20.0%)** | 884 | 9, 43, 164, 199 aa |
| Influenza H3N2 HA | 45 | 3 (6.7%) | 1,195 | 179 aa |

Full listing in `root_length_anomalies.json`.

## They are real prefixes, not misplaced fragments

The truncated roots begin at the canonical Spike start and align cleanly to
position 0 of their leaves:

```
group_005: root 160 aa, leaf 1273 aa, 160/160 root positions mapped
   root[0] -> leaf[0] ; root[159] -> leaf[159]
   root head: MFVFLVLLPLVSSQCVNLTTRTQLPSAYTNSFTRGVYYPDKVFRS
   leaf head: MFVFLVLLPLVSSQCVNLTNRTQLPSAYTNSFTRGVYYPDKVFRS
```

So the alignment fix handles them correctly — nothing is being read at a shifted
offset. The problem is what is *missing*.

## Why it matters more than the frame bug

Generation starts at the root and the models only substitute, so **every
generated leaf inherits the root's length.** Confirmed directly: for
`group_009`, whose root is 160 aa, all 498 generated sequences are 160 aa.

Consequences, in increasing order of severity:

1. **Coverage of the protein.** For these trees only the first ~160 of 1273
   positions are ever scored. The metric silently reports on 12% of Spike.
2. **The RBD is structurally absent.** Spike's receptor-binding domain spans
   residues ~331–531. A 160-aa root cannot contain it, so no RBD mutation can be
   generated, recovered, or scored in 28% of COVID trees. The same trees are the
   ones feeding antigenic recall, the RBD hotspot mask, and the EVEscape tensor —
   every metric the paper leans on for biological relevance.
3. **Training.** 94 of 335 training trees are affected, so this is not only a
   reporting artifact; the bridge is partly fit on truncated Spike.
4. **Coverage and distance metrics too.** `hamming` and `identity` truncate to
   `min(len(a), len(b))`, so a 160-aa generated leaf compared to a 1273-aa
   observed leaf is scored on 160 positions and can look deceptively close.
   This reaches Coverage@K and min-distance — the metrics otherwise untouched by
   the frame fix.

Point 4 means the "SAFE" column in
[`PAPER_NUMBERS_TO_UPDATE.md`](PAPER_NUMBERS_TO_UPDATE.md) is safe *from the
frame bug* but not from this.

## Cause — confirmed

`src/treeencoder/seq_utils.py::nt_to_aa`, called by `run_all_groups.stage_translate`:

```8:21:src/treeencoder/seq_utils.py
    seq = nt_seq.upper().replace('-', '')
    ...
    stop_idx = aa.find('*')
    if stop_idx != -1:
        aa = aa[:stop_idx]
```

augur's ancestral reconstruction emits nucleotides in **alignment coordinates**,
gaps included. Stripping those gaps shifts the reading frame whenever their count
is not a multiple of three; the shifted frame hits a premature stop codon, and
the next line silently truncates the protein there.

This is the same root cause as the frame bug — `.replace('-', '')` applied to a
sequence that is only meaningful in alignment coordinates.

Verified directly on the reconstructed nucleotides:

| tree | gaps in root (nt) | mod 3 | gap-stripped | translated in-frame |
|---|---|---|---|---|
| group_005 | 2 | 2 | **160 aa** | **1273 aa** |
| group_009 | 2 | 2 | **160 aa** | **1273 aa** |
| group_026 | 4 | 1 | **210 aa** | **1273 aa** |

A **two-nucleotide gap destroys 87% of the sequence.** Leaves are unaffected
because real biological indels are multiples of three, so they stay in frame.

## Why this is cheap to fix

**The full-length roots are already on disk.** `group_*_anc_nt.fasta` holds the
correct reconstruction; only the translation stage corrupts it. Fixing this needs
no re-alignment, no FastTree, no `augur refine`, and no `augur ancestral` — just
a corrected re-run of `stage_translate` over files that already exist. Minutes of
CPU.

## What to do

1. **Fix `nt_to_aa`.** Translate in alignment coordinates (codon-wise, gap codons
   to `X`/`-`) instead of stripping gaps first. Stop truncation should only apply
   at a terminal stop; an internal one now means something is wrong and should
   raise rather than silently shorten the protein.
2. **Re-translate**, then assert every root is within a residue or two of the
   alignment width so this cannot recur silently.
3. **Regenerate the affected evaluation roots.** New roots mean new generations,
   so the coverage caches for those roots are invalid. This is GPU work, but only
   for the ~28% of roots affected.
4. **Retrain** to collect the recovered supervision — see below.
5. **Do all of this before the rescoring pass**, or the corrected numbers in
   Tables 3–6 get recomputed on the same broken subset.

## What a retrain buys

The loss masks padding (`valid_mask = (targets != PAD_IDX) & (aa_t != PAD_IDX)`
in `src/bridge/losses.py`), so truncated roots were never supervising the wrong
residues — those positions were simply excluded. The current checkpoint is
**undertrained, not miscalibrated**. Measured over the COVID corpus:

| | supervised positions | RBD positions (331–531) |
|---|---|---|
| train | 77.2M of 102.7M — **24.8% lost** | 11.8M of 16.2M — **27.4% lost** |
| test | 11.2M of 14.9M — 24.6% lost | 1.7M of 2.3M — 27.8% lost |

So a retrain recovers about a quarter of the training signal, and slightly more
than that within the RBD specifically. Past `treesbm_train` jobs ran 3.5–10 h, so
this is affordable.

Worth noting the direction: this is the one correction in the audit expected to
push mutation recall and `AA|hit` **up**, partly offsetting the downward movement
from the frame fix.
