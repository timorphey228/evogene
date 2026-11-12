"""ORF survey for ancestral genomic sequences; complement, both strands, three frames.

ORFs are candidate coding potential only, NOT proof that an ancestor translated them.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict

COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')
STOP = {'TAA', 'TAG', 'TGA'}

@dataclass(frozen=True)
class ORF:
    strand: str
    frame: int
    start: int   # 0-based half-open genomic coordinates
    end: int
    aa_length: int  # excludes stop codon

    def to_dict(self):
        return asdict(self)


def find_orfs(sequence: str, min_aa: int = 30) -> list[ORF]:
    seq = ''.join(sequence.upper().split()).replace('U', 'T')
    if not seq or any(c not in 'ACGTN' for c in seq):
        raise ValueError('Expected DNA letters A/C/G/T/N')
    if min_aa < 1:
        raise ValueError('min_aa must be >=1')
    n = len(seq)
    found: list[ORF] = []
    for strand, s in (('+', seq), ('-', seq.translate(COMPLEMENT)[::-1])):
        for frame in range(3):
            starts: list[int] = []
            for pos in range(frame, n - 2, 3):
                codon = s[pos:pos + 3]
                if codon == 'ATG':
                    starts.append(pos)
                if codon in STOP:
                    for start in starts:
                        length = (pos - start) // 3
                        if length >= min_aa:
                            a, b = (start, pos + 3) if strand == '+' else (n - (pos + 3), n - start)
                            found.append(ORF(strand, frame, a, b, length))
                    starts = []
    return sorted(found, key=lambda x: (-x.aa_length, x.start, x.strand))
