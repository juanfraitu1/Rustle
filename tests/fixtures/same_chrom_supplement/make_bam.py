#!/usr/bin/env python3
"""Generate a small fixture BAM for the same-chromosome supplement integration test.

Design:
- c1 (600 bp) and c2 (320 bp) with block-repeat sequence so canonical GT-AG splice
  sites occur at predictable boundaries.
- Four loci:
    c1:A  start=0,   donor=59,  acceptor=201, end=260  (118 bp spliced, +)
    c2:X  start=0,   donor=59,  acceptor=201, end=260  (118 bp spliced, +)
    c1:B  start=250, donor=299, acceptor=371, end=460  (138 bp spliced, +)
    c1:C  start=300, donor=379, acceptor=531, end=590  (138 bp spliced, +)
- Cross-chrom family: 3 reads (read_cross_1/2/3) with primary on c1:A and secondary
  on c2:X.  Dedicated primary reads on c1:A and c2:X provide skeleton support.
- Same-chrom supplement: 3 reads (read_same_1/2/3) with primary on c1:B and
  secondary on c1:C.  Dedicated primary reads on c1:B and c1:C provide skeleton
  support.
"""
import pysam
import os

# Genome sequence: 10-bp blocks of A, G, T, C repeated.
def genome_seq(n_blocks):
    return "".join(["A" * 10 + "G" * 10 + "T" * 10 + "C" * 10 for _ in range(n_blocks)])

C1_LEN = 600
C2_LEN = 320

ref = {
    "c1": genome_seq(C1_LEN // 40),
    "c2": genome_seq(C2_LEN // 40),
}

# Ensure lengths are exact (n_blocks may overshoot).
ref["c1"] = ref["c1"][:C1_LEN]
ref["c2"] = ref["c2"][:C2_LEN]

fasta_path = os.path.join(os.path.dirname(__file__), "genome.fa")
sam_path = os.path.join(os.path.dirname(__file__), "reads.sam")
bam_path = os.path.join(os.path.dirname(__file__), "reads.bam")

with open(fasta_path, "w") as fh:
    for name, seq in ref.items():
        fh.write(f">{name}\n{seq}\n")

# ---- Locus definitions (0-based) ----
# (chrom, start, donor, acceptor, end, cigar)
LOCI = {
    "c1A": ("c1", 0, 59, 211, 260, "59M152N49M"),
    "c2X": ("c2", 0, 59, 211, 260, "59M152N49M"),
    "c1B": ("c1", 250, 299, 371, 460, "49M72N89M"),
    "c1C": ("c1", 380, 419, 451, 550, "39M32N99M"),
}


def exon_seq(chrom, start, donor, acceptor, end):
    seq = ref[chrom]
    return seq[start:donor] + seq[acceptor:end]


# Pre-compute spliced read sequence for each locus.
LOCI_SEQ = {
    name: exon_seq(chrom, start, donor, acceptor, end)
    for name, (chrom, start, donor, acceptor, end, _) in LOCI.items()
}


def make_record(read_name, chrom, start, cigar, seq, flag, mapq=60, as_score=100, de=0.0, mate_ref="*", mate_pos=0, tlen=0):
    # pysam.AlignedSegment
    a = pysam.AlignedSegment()
    a.query_name = read_name
    a.flag = flag
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigarstring = cigar
    a.next_reference_start = mate_pos
    a.template_length = tlen
    a.query_sequence = seq
    a.query_qualities = None
    a.set_tag("AS", as_score)
    a.set_tag("de", de)
    # reference_id is set later using the output file header.
    return a, chrom, mate_ref


header = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [
        {"SN": "c1", "LN": C1_LEN},
        {"SN": "c2", "LN": C2_LEN},
    ],
}

reads = []

# Dedicated primary reads for each locus (provide skeleton support).
for i in range(3):
    reads.append(make_record(f"c1A_primary_{i}", "c1", 0, LOCI["c1A"][5], LOCI_SEQ["c1A"], 0))
    reads.append(make_record(f"c2X_primary_{i}", "c2", 0, LOCI["c2X"][5], LOCI_SEQ["c2X"], 0))
    reads.append(make_record(f"c1B_primary_{i}", "c1", 250, LOCI["c1B"][5], LOCI_SEQ["c1B"], 0))
    reads.append(make_record(f"c1C_primary_{i}", "c1", 380, LOCI["c1C"][5], LOCI_SEQ["c1C"], 0))

# Cross-chrom reads: primary on c1:A, secondary on c2:X.
for i in range(3):
    reads.append(make_record(f"read_cross_{i}", "c1", 0, LOCI["c1A"][5], LOCI_SEQ["c1A"], 0))
    reads.append(make_record(f"read_cross_{i}", "c2", 0, LOCI["c2X"][5], LOCI_SEQ["c1A"], 256, mapq=0))

# Same-chrom reads: primary on c1:B, secondary on c1:C.
# Use c1:B sequence for the read; the secondary on c1:C will have mismatches but
# the de tag is set to 0.0 so it still ties for the conflict criterion.
for i in range(3):
    reads.append(make_record(f"read_same_{i}", "c1", 250, LOCI["c1B"][5], LOCI_SEQ["c1B"], 0))
    reads.append(make_record(f"read_same_{i}", "c1", 380, LOCI["c1C"][5], LOCI_SEQ["c1B"], 256, mapq=0))

with pysam.AlignmentFile(sam_path, "w", header=header) as outf:
    for rec in reads:
        r, chrom, mate_ref = rec
        r.reference_id = outf.get_tid(chrom)
        if mate_ref != "*":
            r.next_reference_id = outf.get_tid(mate_ref)
        outf.write(r)

pysam.sort("-o", bam_path, sam_path)
pysam.index(bam_path)
# Regenerate FASTA index now that the FASTA has been (re)written.
import subprocess
subprocess.run(["samtools", "faidx", fasta_path], check=True)
print(f"Wrote {bam_path}")
