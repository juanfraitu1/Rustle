#!/usr/bin/env python3
"""Generate the `o2_origin_resolve` regression fixture.

The fixture proves that copy_assign's O2 PSV-resolution overrides the
aligner's primary flag: a molecule whose primary record sits at copy A but
whose sequence matches copy B is assigned to copy B, the assembled transcript
is built at copy B, and `missing_copy_flag --assignments` attributes the read
to copy B.

Layout (one contig `c1`, 5000 bp; both copies are two-exon, 600 bp spliced
length, `+`):

    copy A  genomic c1:500-1200   spliced cDNA pattern 0
    copy B  genomic c1:2600-3300  spliced cDNA pattern 1

Each copy has exons 100 bp + 500 bp separated by a 100 bp canonical intron.
The two spliced cDNA patterns differ at 24 PSV positions over the 600 bp
spliced length, more than enough for O2 to distinguish them.

Support reads at each copy start at four different 5' offsets inside the
first exon but share the same intron chain, so they survive the
coordinate-based primary de-duplication and are collapsed into one isoform by
the assembler.

Planted molecules:
  * 4 support reads at copy A carrying pattern 0.
  * 4 support reads at copy B carrying pattern 1.
  * `MOL_PRIMARY_WRONG` -- one molecule with two records:
      - primary   at copy A (offset 5) carrying pattern 1 (mismatches copy A).
      - secondary at copy B (offset 5) carrying pattern 1 (matches copy B).
      - both records have the same AS score (= aligned query length) so the
        default AS-tied gate admits the molecule.
"""
import os
import subprocess

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
CHROM = "c1"
CHROM_LEN = 5000
COPY_GENOMIC_LEN = 700
E1_LEN = 100
E2_LEN = 500
INTRON_LEN = 100
N_SUPPORT = 4
BASES = "ACGT"

# 24 planted PSV offsets inside the 600 bp spliced cDNA.
PSV_OFF = [20 + 24 * j for j in range(24)]


def lcg(seed):
    """Deterministic pseudo-random byte stream."""
    x = seed & 0xFFFFFFFF
    while True:
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        yield (x >> 16) & 0xFF


def filler(n, seed):
    g = lcg(seed)
    return "".join(BASES[next(g) & 3] for _ in range(n))


CDNA_LEN = E1_LEN + E2_LEN
BACKBONE = filler(CDNA_LEN, 7)


def pattern(k):
    """Copy pattern k: shared backbone with a k-specific base at every PSV."""
    s = list(BACKBONE)
    g = lcg(1000 + 37 * k)
    for off in PSV_OFF:
        s[off] = BASES[next(g) & 3]
    return "".join(s)


PAT = [pattern(k) for k in range(2)]
for a in range(2):
    for b in range(a + 1, 2):
        d = sum(1 for i in range(CDNA_LEN) if PAT[a][i] != PAT[b][i])
        assert d >= 8, f"patterns {a},{b} differ at only {d} positions"


def genomic_copy(seq_list, start, pat):
    """Embed a two-exon copy into the genome list (modifies in place)."""
    seq_list[start : start + E1_LEN] = list(pat[:E1_LEN])
    intron = [BASES[(i + start) & 3] for i in range(INTRON_LEN)]
    intron[0] = "G"
    intron[1] = "T"
    intron[-2] = "A"
    intron[-1] = "G"
    seq_list[start + E1_LEN : start + E1_LEN + INTRON_LEN] = intron
    seq_list[start + E1_LEN + INTRON_LEN : start + COPY_GENOMIC_LEN] = list(
        pat[E1_LEN:]
    )


# ---- genome ---------------------------------------------------------------------------------
A_START = 500
B_START = 2600
seq = list(filler(CHROM_LEN, 99))
genomic_copy(seq, A_START, PAT[0])
genomic_copy(seq, B_START, PAT[1])
REF = "".join(seq)

with open(os.path.join(HERE, "genome.fa"), "w") as fh:
    fh.write(f">{CHROM}\n")
    for i in range(0, len(REF), 60):
        fh.write(REF[i : i + 60] + "\n")

# ---- catalog --------------------------------------------------------------------------------
with open(os.path.join(HERE, "copies.tsv"), "w") as fh:
    fh.write(
        "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n"
    )
    for ci, start in enumerate([A_START, B_START]):
        end = start + COPY_GENOMIC_LEN
        e1_end = start + E1_LEN
        e2_start = start + E1_LEN + INTRON_LEN
        fh.write(
            f"FAM1\t{ci}\tDN_{CHROM}_{start}_2\t{CHROM}\t{start}\t{end}\t2\t+\t{N_SUPPORT}\t"
            f"{start}-{e1_end},{e2_start}-{end}\n"
        )

with open(os.path.join(HERE, "copies.fa"), "w") as fh:
    for ci, start in enumerate([A_START, B_START]):
        end = start + COPY_GENOMIC_LEN
        fh.write(f">FAM1|{ci}|{CHROM}:{start}-{end}|+|nexon=2\n{PAT[ci]}\n")

# ---- BAM ------------------------------------------------------------------------------------
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
records = []


def spliced_read(copy_start, offset, total_cdna_len):
    """Return a spliced read aligned to a copy.

    The read starts `offset` bases into the spliced cDNA and covers
    `total_cdna_len` bases of the cDNA.  It is aligned as one block in exon 1,
    one intron-sized N, and one block in exon 2.  All reads from the same copy
    share the same intron chain so the assembler groups them, while their
    different 5' offsets keep them distinct through the coordinate de-duplication
    of `primary`.
    """
    genomic_start = copy_start + offset
    e1_aligned = E1_LEN - offset
    e2_aligned = total_cdna_len - e1_aligned
    ref_end = copy_start + E1_LEN + INTRON_LEN + e2_aligned
    cigar = f"{e1_aligned}M{INTRON_LEN}N{e2_aligned}M"
    return genomic_start, ref_end, cigar


def rec(name, copy_start, offset, total_cdna_len, seq, flag, mapq=60):
    genomic_start, ref_end, cigar = spliced_read(copy_start, offset, total_cdna_len)
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    a.reference_start = genomic_start
    a.mapping_quality = mapq
    a.cigarstring = cigar
    a.query_sequence = seq[offset : offset + total_cdna_len]
    a.query_qualities = None
    a.set_tag("AS", total_cdna_len)
    a.set_tag("de", 0.0)
    records.append(a)


# support reads: four different 5' offsets inside exon 1, all running to the
# end of the spliced cDNA.  They share the intron chain of their copy.
SUPPORT_OFFSETS = [0, 15, 30, 45]
for i, off in enumerate(SUPPORT_OFFSETS):
    rec(f"supA_{i}", A_START, off, CDNA_LEN - off, PAT[0], 0)
    rec(f"supB_{i}", B_START, off, CDNA_LEN - off, PAT[1], 0)

# the molecule whose primary is at copy A but whose sequence belongs to copy B
WRONG_OFF = 5
rec("MOL_PRIMARY_WRONG", A_START, WRONG_OFF, CDNA_LEN - WRONG_OFF, PAT[1], 0)
rec("MOL_PRIMARY_WRONG", B_START, WRONG_OFF, CDNA_LEN - WRONG_OFF, PAT[1], 256, mapq=0)

sam_path = os.path.join(HERE, "reads.sam")
bam_path = os.path.join(HERE, "reads.bam")
with pysam.AlignmentFile(sam_path, "w", header=header) as out:
    for a in records:
        a.reference_id = out.get_tid(CHROM)
        out.write(a)

pysam.sort("-o", bam_path, sam_path)
pysam.index(bam_path)
subprocess.run(["samtools", "faidx", os.path.join(HERE, "genome.fa")], check=True)
print(f"wrote {bam_path} ({len(records)} records)")
