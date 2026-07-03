"""Nucleotide → amino acid translation utilities for influenza HA sequences."""

from Bio.Seq import Seq


def nt_to_aa(nt_seq: str, cds_start: int | None = None) -> str:
    """
    Translate a nucleotide sequence to amino acids.

    Finds the first ATG automatically if cds_start is not given.
    Returns the AA sequence up to (not including) the first stop codon.

    Args:
        nt_seq:    nucleotide string (may contain leading/trailing UTR)
        cds_start: position of the first ATG; auto-detected if None

    Returns:
        amino acid string (no stop codon)
    """
    seq = nt_seq.upper().replace('-', '')

    if cds_start is None:
        cds_start = seq.find('ATG')
        if cds_start == -1:
            raise ValueError("No ATG start codon found in sequence")

    cds = seq[cds_start:]
    # Trim to multiple of 3
    cds = cds[:len(cds) - len(cds) % 3]
    aa = str(Seq(cds).translate())

    # Trim at first stop codon
    stop_idx = aa.find('*')
    if stop_idx != -1:
        aa = aa[:stop_idx]

    return aa


def translate_asr_fasta(input_fasta: str, output_fasta: str) -> int:
    """
    Translate all sequences in an ASR nucleotide FASTA to amino acids.

    Returns number of sequences written.
    """
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    records = []
    for rec in SeqIO.parse(input_fasta, 'fasta'):
        aa = nt_to_aa(str(rec.seq))
        records.append(SeqRecord(Seq(aa), id=rec.id, description=''))

    SeqIO.write(records, output_fasta, 'fasta')
    return len(records)
