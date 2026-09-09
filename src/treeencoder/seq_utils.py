#Nucleotide to amino acid translation utilities for influenza HA sequences.

from Bio.Seq import Seq

# augur's ancestral reconstruction emits nucleotides in *alignment* coordinates,
# gaps included. Stripping those gaps before translating shifts the reading frame
# whenever their count is not a multiple of three; the shifted frame then hits a
# premature stop codon and the protein gets truncated there. On SARS-CoV-2 a
# two-nucleotide gap in the reconstructed root truncated Spike from 1273 aa to
# 160 aa, and it did so in 28% of trees.
#
# The inputs here are in-frame CDS slices, so three alignment columns are one
# reference codon and translating column-wise keeps the frame. See
# results/voc_threat_panel/TRUNCATED_ROOTS.md.


def _frameshifting_gaps(seq: str) -> bool:
    """True if stripping this sequence's gaps would shift the reading frame.

    A real biological indel is a whole number of codons, so every maximal gap run
    has a length divisible by three and removing them preserves the frame. A run
    that is not divisible by three cannot be a coding indel; it is an alignment
    or reconstruction artifact, and stripping it corrupts everything downstream.
    """
    run = 0
    for ch in seq:
        if ch == '-':
            run += 1
        else:
            if run % 3:
                return True
            run = 0
    return run % 3 != 0


def _find_cds_start_aligned(seq: str) -> int:
    """Alignment-coordinate index of the first ATG."""
    ungapped = seq.replace('-', '')
    i = ungapped.find('ATG')
    if i == -1:
        raise ValueError("No ATG start codon found in sequence")
    seen = 0
    for k, ch in enumerate(seq):
        if ch == '-':
            continue
        if seen == i:
            return k
        seen += 1
    raise ValueError("No ATG start codon found in sequence")


def _drop_terminal_stop(aa: str) -> str:
    """Remove a trailing stop; convert any internal stop to ``X``.

    A stop followed only by gaps or ambiguity is the real end of the CDS. A stop
    with actual protein after it means the frame or the reconstruction is wrong,
    and truncating there is what shortened HIV roots to a fifth of their leaves'
    length. Keeping the residue count intact and flagging the position as
    unknown loses far less than deleting 80% of the protein.
    """
    stop_idx = aa.find('*')
    if stop_idx == -1:
        return aa
    trailing = aa[stop_idx + 1:].replace('-', '').replace('X', '')
    if not trailing:
        return aa[:stop_idx]
    return aa.replace('*', 'X')


def nt_to_aa(
    nt_seq: str,
    cds_start: int | None = None,
    mode: str = "aligned",
) -> str:
    """Translate a (possibly aligned) nucleotide CDS to amino acids.

    ``mode`` selects how alignment gaps are treated, and the right choice is a
    property of the *alignment*, not of the individual sequence:

    ``"aligned"``
        Three alignment columns are one reference codon, so translating
        column-wise preserves the frame: an all-gap codon is a real deletion and
        is dropped, a partial-gap codon is the frameshift artifact and becomes
        ``X``. Correct only for reference-anchored alignments. Sequences whose
        gap runs are all codon-sized take the plain path below, so their output
        is byte-identical to the pre-fix translator.

    ``"ungapped"``
        Gaps carry no cross-sequence frame, so each sequence is stripped and
        translated in its own frame. This is right for per-group MAFFT
        nucleotide alignments, where no single column offset is in frame for
        every sequence and column-wise translation would manufacture ``X`` runs.

    ``"legacy"``
        Strip gaps and truncate at the first stop. Reproduces pre-fix output;
        it is the bug, kept only for regression comparisons.

    See scripts/panviral/audit_frame_offset.py for the evidence behind each
    dataset's assignment.
    """
    if mode not in ("aligned", "ungapped", "legacy"):
        raise ValueError(f"unknown mode {mode!r}")

    seq = nt_seq.upper()

    if mode != "aligned" or not _frameshifting_gaps(seq):
        seq = seq.replace('-', '')
        if cds_start is None:
            cds_start = seq.find('ATG')
            if cds_start == -1:
                raise ValueError("No ATG start codon found in sequence")
        cds = seq[cds_start:]
        cds = cds[:len(cds) - len(cds) % 3]
        aa = str(Seq(cds).translate())
        if mode == "legacy":
            stop_idx = aa.find('*')
            return aa[:stop_idx] if stop_idx != -1 else aa
        return _drop_terminal_stop(aa)

    if cds_start is None:
        cds_start = _find_cds_start_aligned(seq)

    out = []
    for j in range(cds_start, len(seq) - 2, 3):
        codon = seq[j:j + 3]
        n_gap = codon.count('-')
        if n_gap == 3:
            out.append('-')
        elif n_gap:
            out.append('X')
        else:
            try:
                out.append(str(Seq(codon).translate()))
            except Exception:  # noqa: BLE001 - ambiguous nucleotides
                out.append('X')

    aa = ''.join(out)

    # A terminal stop is expected; drop it (and anything after it, which is
    # 3'UTR read-through rather than protein).
    stop_idx = aa.find('*')
    if stop_idx != -1:
        trailing = aa[stop_idx + 1:].replace('-', '').replace('X', '')
        if not trailing:
            aa = aa[:stop_idx]
        else:
            # An internal stop in a correct frame is a genuine anomaly. Keep the
            # protein full length rather than silently truncating it -- silent
            # truncation is the bug this function exists to avoid.
            aa = aa.replace('*', 'X')

    return aa.replace('-', '')


def translate_asr_fasta(input_fasta: str, output_fasta: str, gap_aware: bool = True) -> int:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    records = []
    for rec in SeqIO.parse(input_fasta, 'fasta'):
        aa = nt_to_aa(str(rec.seq), gap_aware=gap_aware)
        records.append(SeqRecord(Seq(aa), id=rec.id, description=''))

    SeqIO.write(records, output_fasta, 'fasta')
    return len(records)
