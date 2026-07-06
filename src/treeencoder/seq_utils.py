#Nucleotide to amino acid translation utilities for influenza HA sequences.

from Bio.Seq import Seq


def nt_to_aa(nt_seq: str, cds_start: int | None = None) -> str:
    #translates a nucleotide sequence to amino acids, starting at the first ATG if cds_start is None
    seq = nt_seq.upper().replace('-', '')

    if cds_start is None:
        cds_start = seq.find('ATG')
        if cds_start == -1:
            raise ValueError("No ATG start codon found in sequence")

    cds = seq[cds_start:]
    cds = cds[:len(cds) - len(cds) % 3]
    aa = str(Seq(cds).translate())

    stop_idx = aa.find('*')
    if stop_idx != -1:
        aa = aa[:stop_idx]

    return aa


def translate_asr_fasta(input_fasta: str, output_fasta: str) -> int:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    records = []
    for rec in SeqIO.parse(input_fasta, 'fasta'):
        aa = nt_to_aa(str(rec.seq))
        records.append(SeqRecord(Seq(aa), id=rec.id, description=''))

    SeqIO.write(records, output_fasta, 'fasta')
    return len(records)
