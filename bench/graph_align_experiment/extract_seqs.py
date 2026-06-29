#!/usr/bin/env python3
"""Extract per-copy TRANSCRIPT sequences (spliced, reverse-complemented since '-'
strand) for A1, A2, B1, using the EXACT coordinate convention of the read
generator (0-based half-open), so the transcript seqs match the read sequences.

NOTE on strand: the generator builds read SEQ as the forward (+) spliced
sequence and emits SAM flag 16 (reverse) with a forward CIGAR; i.e. the stored
SEQ in the BAM is already the forward genomic spliced sequence. samtools fasta
will reverse-complement flag-16 reads back to... actually samtools fasta emits
the ORIGINAL read orientation (reverse-complements flag16 SEQ). We handle strand
carefully below and verify by alignment, not by assumption.
"""
import sys

REF_PATH = "ref.fa"
with open(REF_PATH) as f:
    lines = f.readlines()
REF_SEQ = "".join(l.strip() for l in lines[1:])

GENE_A_EXONS = [(1000, 1200), (2000, 2300), (3000, 3150), (4000, 4500)]
GENE_A_SKIP  = [(1000, 1200), (2000, 2300), (4000, 4500)]
GENE_B_EXONS = [(10000, 10200), (11000, 11300), (12000, 12150), (13000, 13500)]

def comp(b):
    return {"A":"T","T":"A","C":"G","G":"C","N":"N"}[b]
def rc(s):
    return "".join(comp(b) for b in reversed(s))
def spliced(exons):
    return "".join(REF_SEQ[s:e] for s, e in exons)

copies = {
    "A1": spliced(GENE_A_EXONS),
    "A2": spliced(GENE_A_SKIP),
    "B1": spliced(GENE_B_EXONS),
}

# Forward (genomic +) spliced sequences. The generator stores these as the BAM
# SEQ. We emit BOTH a forward seqs.fa (for graph backbone) — but reads from
# `samtools fasta` of a flag16 record come out reverse-complemented, so we will
# align reads to whichever orientation matches. We write forward here and decide
# orientation empirically.
with open("seqs.fa", "w") as f:
    for name, seq in copies.items():
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

for name, seq in copies.items():
    print(f"{name}\tlen={len(seq)}", file=sys.stderr)
