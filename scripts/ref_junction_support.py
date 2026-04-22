"""
Index splice junctions from a BAM (primary alignments) and match reference GFF introns.

Coordinate convention matches src/bam_parser.rs (CIGAR N / long D): donor and acceptor
are 1-based genomic positions — last base of upstream exon and first base of downstream exon.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterator, List, Tuple

try:
    import pysam
except ImportError as e:
    raise SystemExit("ref_junction_support requires pysam") from e

# CIGAR operation codes (pysam)
CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_EQ = 7
CIGAR_X = 8


def read_strand(read: pysam.AlignedSegment) -> str:
    xs = read.get_tag("XS") if read.has_tag("XS") else None
    if xs in (b"+", "+"):
        return "+"
    if xs in (b"-", "-"):
        return "-"
    return "-" if read.is_reverse else "+"


def iter_read_junctions_1bp(
    read: pysam.AlignedSegment,
) -> Iterator[Tuple[str, str, int, int]]:
    """Yield (chrom, strand, donor_1b, acceptor_1b) for each splice in a primary alignment."""
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return
    chrom = read.reference_name
    if chrom is None:
        return
    strand = read_strand(read)
    t = read.cigartuples
    if not t:
        return
    ref_pos = read.reference_start
    for op, ln in t:
        if op in (CIGAR_M, CIGAR_EQ, CIGAR_X):
            ref_pos += ln
        elif op == CIGAR_D:
            if 5 <= ln <= 1_000_000:
                donor = ref_pos
                ref_pos += ln
                acceptor = ref_pos + 1
                yield chrom, strand, donor, acceptor
            else:
                ref_pos += ln
        elif op == CIGAR_N:
            donor = ref_pos
            ref_pos += ln
            acceptor = ref_pos + 1
            yield chrom, strand, donor, acceptor
        elif op in (CIGAR_I, CIGAR_S, CIGAR_H):
            continue


def quantize_junction(d: int, a: int, slop: int) -> Tuple[int, int]:
    if slop <= 0:
        return d, a
    w = 2 * slop + 1
    return d // w, a // w


def junction_store_key(
    chrom: str, strand: str, d: int, a: int, slop: int
) -> Tuple[str, str, int, int]:
    dq, aq = quantize_junction(d, a, slop)
    return chrom, strand, dq, aq


def build_junction_support_counter(
    bam: pysam.AlignmentFile, slop: int
) -> Counter[Tuple[str, str, int, int]]:
    """Count primary spliced junctions (optionally quantized for tolerance)."""
    ctr: Counter[Tuple[str, str, int, int]] = Counter()
    for read in bam:
        for chrom, strand, d, a in iter_read_junctions_1bp(read):
            ctr[junction_store_key(chrom, strand, d, a, slop)] += 1
    return ctr


def build_junction_support_counter_from_path(
    bam_path, slop: int
) -> Counter[Tuple[str, str, int, int]]:
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return build_junction_support_counter(bam, slop)


@dataclass
class RefJunctionStats:
    num_introns: int
    introns_with_support: int
    min_junction_reads: int  # 0 if num_introns == 0; else min over introns
    full_chain_supported: bool  # True if num_introns == 0 or all supported


def ref_introns_from_exons(
    exons: List[Tuple[int, int]], strand: str
) -> List[Tuple[int, int]]:
    """Sorted by genomic start; return (donor_1b, acceptor_1b) per intron."""
    if len(exons) < 2:
        return []
    se = sorted(exons)
    out: List[Tuple[int, int]] = []
    for i in range(len(se) - 1):
        _s1, e1 = se[i]
        s2, _e2 = se[i + 1]
        if s2 <= e1:
            continue
        out.append((e1, s2))
    return out


def count_junction_support_for_transcript(
    counter: Counter[Tuple[str, str, int, int]],
    chrom: str,
    strand: str,
    introns: List[Tuple[int, int]],
    slop: int,
) -> RefJunctionStats:
    if not introns:
        return RefJunctionStats(0, 0, 0, True)
    counts: List[int] = []
    supported = 0
    for d, a in introns:
        k = junction_store_key(chrom, strand, d, a, slop)
        c = counter.get(k, 0)
        counts.append(c)
        if c > 0:
            supported += 1
    mmin = min(counts)
    return RefJunctionStats(
        len(introns), supported, mmin, supported == len(introns)
    )


def parse_gff_exons_by_transcript(
    path,
) -> Tuple[Dict[str, Tuple[str, str, List[Tuple[int, int]]]], Dict[str, int]]:
    """transcript_id -> (chrom, strand, exons). Also num_exons per tid (for bins)."""
    import re
    from collections import defaultdict

    exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    meta: Dict[str, Tuple[str, str]] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            if p[2] not in ("transcript", "exon"):
                continue
            m = re.search(r'transcript_id\s+"([^"]+)"', p[8])
            if not m:
                continue
            tid = m.group(1)
            if p[2] == "transcript":
                meta[tid] = (p[0], p[6])
            else:
                exons[tid].append((int(p[3]), int(p[4])))
    out: Dict[str, Tuple[str, str, List[Tuple[int, int]]]] = {}
    nexon: Dict[str, int] = {}
    for tid, (chrom, strand) in meta.items():
        ev = exons.get(tid, [])
        if not ev:
            continue
        out[tid] = (chrom, strand, ev)
        nexon[tid] = len(ev)
    return out, nexon
