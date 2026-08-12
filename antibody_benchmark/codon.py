"""Minimal NT->AA helpers for shared metrics (standard genetic code)."""

from __future__ import annotations

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_nt(seq: str, stop_as: str = "X", ambig_as: str = "X") -> str:
    """Translate nucleotide sequence; non-ACGT codons -> ambig_as; stops -> stop_as."""
    s = seq.upper().replace("U", "T")
    out = []
    for i in range(0, len(s) - len(s) % 3, 3):
        codon = s[i : i + 3]
        if any(b not in "ACGT" for b in codon):
            out.append(ambig_as)
            continue
        aa = CODON_TABLE.get(codon, ambig_as)
        out.append(stop_as if aa == "*" else aa)
    return "".join(out)


def looks_like_nt(seq: str) -> bool:
    chars = set(seq.upper()) - {"-", "."}
    return bool(chars) and chars <= set("ACGTNXY")
