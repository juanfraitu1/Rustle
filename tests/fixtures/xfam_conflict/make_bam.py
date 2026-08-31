#!/usr/bin/env python3
"""Generate the TWO-FAMILY cross-family conflict fixture (`tests/copy_assign_xfam.rs`).

Why it exists: every pre-existing `copy_assign` e2e fixture supplies exactly ONE family
(`tests/copy_assign_families.rs`'s GWFAM1), so no in-tree test could ever observe a molecule
being ACCEPTED by two families' significance gates at once — which is precisely why the
cross-family double-assignment defect shipped. This fixture supplies TWO families on one contig
and plants one molecule per stratum of the reconciliation rule.

Layout (contig `x1`, 20 kb; all copies single-exon, 600 bp, `+`):

    XFA copy0  x1:1000-1600   pattern P0
    XFA copy1  x1:3000-3600   pattern P1
    XFA copy2  x1:15000-15600 pattern P2      <-- same interval as XFB copy2
    XFB copy0  x1:8000-8600   pattern P4
    XFB copy1  x1:10000-10600 pattern P3
    XFB copy2  x1:15000-15600 pattern P2      <-- same interval as XFA copy2

The five loci carry five pairwise-distinct PSV patterns over 24 planted variant positions, so
within each family a read carrying one copy's pattern is decisively assignable to it.

Planted molecules (one per stratum):

  MOL_CONTRA       TWO records, one primary at XFA copy0 (x1:1000-1600) carrying P0 and one
                   SECONDARY at XFB copy1 (x1:10000-10600) carrying P3. The two claimed copy
                   intervals are DISJOINT and the two records are DIFFERENT placements, so a
                   single molecule is claimed by two loci it cannot both have come from
                   => cross_family_contradiction (the only stratum that abstains).
  MOL_READTHROUGH  ONE record, x1:3000 with CIGAR 600M4400N600M, sequence P1+P4: its N-gap spans
                   XFA copy1 to XFB copy0. Disjoint copies but ONE record => readthrough_span:
                   the molecule really does span both loci, so the contradiction argument does
                   not apply and it must NOT abstain.
  MOL_SHARED       ONE record inside x1:15000-15600, carrying P2 — the interval XFA copy2 and
                   XFB copy2 BOTH claim => shared_locus: an O1 partition artifact (one locus in
                   two families), not an O2 assignment error, so it must NOT abstain either.

Support reads: 8 primary reads per locus (>= PSV_MIN_ALLELE_READS=2 per allele and
>= PSV_MIN_JUDGE_COV=4 judgeable coverage), and the `--families` contract requires every
supplied copy to have at least one overlapping read.
"""
import os
import subprocess

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
CHROM = "x1"
CHROM_LEN = 20000
COPY_LEN = 600
N_SUPPORT = 8
BASES = "ACGT"

# 24 planted variant offsets inside a 600 bp copy.
PSV_OFF = [20 + 24 * j for j in range(24)]


def lcg(seed):
    """Deterministic pseudo-random byte stream (no numpy / no random module state)."""
    x = seed & 0xFFFFFFFF
    while True:
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        yield (x >> 16) & 0xFF


def filler(n, seed):
    g = lcg(seed)
    return "".join(BASES[next(g) & 3] for _ in range(n))


# ---- five pairwise-distinct copy patterns over one 600 bp backbone ---------------------------
BACKBONE = filler(COPY_LEN, 7)


def pattern(k):
    """Copy pattern k: the shared backbone with a k-specific base at every PSV offset."""
    s = list(BACKBONE)
    g = lcg(1000 + 37 * k)
    for off in PSV_OFF:
        s[off] = BASES[next(g) & 3]
    return "".join(s)


PAT = [pattern(k) for k in range(5)]
# The fixture is only meaningful if the patterns are mutually distinguishable at enough columns.
for a in range(5):
    for b in range(a + 1, 5):
        d = sum(1 for i in range(COPY_LEN) if PAT[a][i] != PAT[b][i])
        assert d >= 8, f"patterns {a},{b} differ at only {d} positions"

# ---- genome ---------------------------------------------------------------------------------
# (start, pattern index) for each planted locus.
LOCI = [(1000, 0), (3000, 1), (8000, 4), (10000, 3), (15000, 2)]
seq = list(filler(CHROM_LEN, 99))
for start, pi in LOCI:
    seq[start : start + COPY_LEN] = list(PAT[pi])
REF = "".join(seq)

with open(os.path.join(HERE, "genome.fa"), "w") as fh:
    fh.write(f">{CHROM}\n")
    for i in range(0, len(REF), 60):
        fh.write(REF[i : i + 60] + "\n")

# ---- catalog (copies.tsv + copies.fa), the `--families` input --------------------------------
# copy_idx follows ascending start, matching gw_family_catalog's own ordering.
CATALOG = [
    ("XFA", 0, 1000, 0),
    ("XFA", 1, 3000, 1),
    ("XFA", 2, 15000, 2),
    ("XFB", 0, 8000, 4),
    ("XFB", 1, 10000, 3),
    ("XFB", 2, 15000, 2),
]

with open(os.path.join(HERE, "copies.tsv"), "w") as fh:
    fh.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
    for fam, ci, start, _pi in CATALOG:
        end = start + COPY_LEN
        fh.write(
            f"{fam}\t{ci}\tDN_{CHROM}_{start}_1_{fam}\t{CHROM}\t{start}\t{end}\t1\t+\t{N_SUPPORT}\t{start}-{end}\n"
        )

with open(os.path.join(HERE, "copies.fa"), "w") as fh:
    for fam, ci, start, pi in CATALOG:
        end = start + COPY_LEN
        fh.write(f">{fam}|{ci}|{CHROM}:{start}-{end}|+|nexon=1\n{PAT[pi]}\n")

# ---- BAM ------------------------------------------------------------------------------------
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
records = []


def rec(name, start, cigar, s, flag, mapq=60, as_score=600):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigarstring = cigar
    a.query_sequence = s
    a.query_qualities = None
    a.set_tag("AS", as_score)
    a.set_tag("de", 0.0)
    records.append(a)


# support reads: N_SUPPORT primaries per planted locus.
for start, pi in LOCI:
    for i in range(N_SUPPORT):
        rec(f"sup_{start}_{i}", start, f"{COPY_LEN}M", PAT[pi], 0)

# stratum 3 — two independent placements naming disjoint loci.
rec("MOL_CONTRA", 1000, f"{COPY_LEN}M", PAT[0], 0)
rec("MOL_CONTRA", 10000, f"{COPY_LEN}M", PAT[3], 256, mapq=0)

# stratum 2 — one record whose N-gap spans XFA copy1 -> XFB copy0.
rec("MOL_READTHROUGH", 3000, f"{COPY_LEN}M4400N{COPY_LEN}M", PAT[1] + PAT[4], 0, as_score=1200)

# stratum 1 — one record inside the interval both families claim.
rec("MOL_SHARED", 15000, f"{COPY_LEN}M", PAT[2], 0)

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
