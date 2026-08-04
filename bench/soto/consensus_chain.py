#!/usr/bin/env python3
"""Build a locus representative from the JUNCTION SET instead of picking one observed intron chain.

Usage: consensus_chain.py <copies.tsv> <bam> [chr1,chr15]

WHY THE EXISTING FIXES ARE MARGINAL. `pick_locus_rep` scores candidates by read support, and unspliced
reads pool into ONE cluster by span overlap while spliced reads are grouped by EXACT intron chain and so
shatter. `RUSTLE_SPLICED_REP` fixed the comparison by summing the spliced chains into one class -- but the
representative is then still a single observed chain, and the largest chain at NOTCH2 carries 5% of the
reads. There is no good chain to pick. Measured on the same loci: NOTCH2NLC has 117 distinct junction sets
from 156 spliced reads.

WHY JUNCTIONS SURVIVE WHAT CHAINS DO NOT. Joining a chain requires a read to match EVERY junction, so chain
support decays combinatorially. A junction is supported independently. At NOTCH2NLC, 55 junctions carry >= 3
reads and 32 carry >= 10, and 29 of 37 stub-represented truth loci have at least one junction at >= 3 reads.

THE CONSTRUCTION. Take every junction at the locus with >= MIN_READS support and choose the maximum-weight
subset that is mutually COMPATIBLE (introns pairwise non-overlapping, so the result is a valid transcript
rather than a chimera of mutually exclusive isoforms). Sorting introns by end coordinate makes this a
longest-path problem on a DAG, solved exactly by one O(n log n) dynamic program. No tie-break, no heuristic
ranking, one threshold.

WHAT IS MEASURED. For each truth member whose current representative is a single-exon stub: does a
consensus chain exist, how many junctions does it carry, and how does its span compare to the true gene
span (which the shipped stub gets wrong by 0.03-0.52).
"""
import subprocess
import sys
from bisect import bisect_right
from collections import defaultdict

COP, BAM = sys.argv[1], sys.argv[2]
CHROMS = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else None
BED = "bench/soto/80_fams.gene_preferred.bed"
MIN_READS = 3

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 7 or (CHROMS and f[3] not in CHROMS):
        continue
    copies.append((f[3], int(f[4]), int(f[5]), int(f[6])))

members = []
for ln in open(BED):
    c, s, e, nm, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and c not in CHROMS:
        continue
    members.append((c, int(s), int(e), nm.split("|")[0]))


def junctions(chrom, lo, hi):
    out = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{chrom}:{lo}-{hi}"],
                         capture_output=True, text=True).stdout
    J = defaultdict(int)
    for ln in out.splitlines():
        f = ln.split("\t")
        p, n = int(f[3]), 0
        for ch in f[5]:
            if ch.isdigit():
                n = n * 10 + ord(ch) - 48
            else:
                if ch == "N":
                    J[(p, p + n)] += 1
                    p += n
                elif ch in "MDX=":
                    p += n
                n = 0
    return J


def max_weight_chain(J, min_reads):
    """Max-weight set of pairwise NON-OVERLAPPING introns. Exact: DP over introns sorted by end."""
    iv = sorted(((a, b, w) for (a, b), w in J.items() if w >= min_reads), key=lambda x: x[1])
    if not iv:
        return []
    ends = [b for _, b, _ in iv]
    best = [0] * (len(iv) + 1)
    take = [None] * len(iv)
    for i, (a, b, w) in enumerate(iv):
        j = bisect_right(ends, a) - 1          # last intron ending at or before this one starts
        with_i = w + (best[j + 1] if j >= 0 else 0)
        if with_i > best[i]:
            best[i + 1], take[i] = with_i, j
        else:
            best[i + 1], take[i] = best[i], None
    chain, i = [], len(iv) - 1
    while i >= 0:
        if best[i + 1] != best[i]:
            chain.append(iv[i])
            i = take[i] if take[i] is not None else -1
        else:
            i -= 1
    return sorted(chain)


rows = []
for (c, ms, me, g) in members:
    ov = [(min(me, pe) - max(ms, ps), ps, pe, nx) for (pc, ps, pe, nx) in copies
          if pc == c and ms < pe and me > ps]
    if not ov:
        continue
    _, ps, pe, nx = max(ov)
    if nx > 1:
        continue                                   # only the stub-represented loci
    J = junctions(c, ps, pe)
    ch = max_weight_chain(J, MIN_READS)
    T = me - ms
    if ch:
        span = ch[-1][1] - ch[0][0]
        rows.append((g, (pe - ps) / T, len(ch), span / T, min(w for _, _, w in ch)))
    else:
        rows.append((g, (pe - ps) / T, 0, 0.0, 0))

n = len(rows)
have = [r for r in rows if r[2] > 0]
print(f"stub-represented truth members: {n}")
print(f"  a consensus chain exists (>= {MIN_READS} reads/junction): {len(have)} "
      f"({100*len(have)/max(n,1):.0f}%)\n")
if have:
    j = sorted(r[2] for r in have)
    sp = sorted(r[3] for r in have)
    st = sorted(r[1] for r in have)
    m = len(have)
    print(f"  junctions in the chain : median {j[m//2]}  max {j[-1]}")
    print(f"  shipped stub span/truth: median {st[m//2]:.2f}   in 0.5-2x "
          f"{100*sum(1 for x in st if 0.5<=x<=2)/m:.0f}%")
    print(f"  chain span / truth     : median {sp[m//2]:.2f}   in 0.5-2x "
          f"{100*sum(1 for x in sp if 0.5<=x<=2)/m:.0f}%")

print(f"\n{'gene':<14}{'stub':>7}{'junc':>6}{'chain':>8}{'minRd':>7}")
for g, st, nj, sp, mw in sorted(rows, key=lambda r: -r[2])[:18]:
    print(f"{g:<14}{st:>7.2f}{nj:>6}{sp:>8.2f}{mw:>7}")
