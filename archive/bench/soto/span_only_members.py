#!/usr/bin/env python3
"""How often does a transcript join a locus by SPAN alone, contributing no exonic sequence to it?

Usage: span_only_members.py <bam> <copies.tsv> <chrom> [start] [end]

Locus membership is decided by span overlap -- `ref_start < rep.end && read_end > rep.start`. A giant
intron satisfies that vacuously, so a transcript that splices straight OVER a gene is counted among its
reads. An intron is precisely the assertion that the transcript is NOT present there, so this is membership
conferred by absence.

Found at NPIPB9, where the 34-read chain that `pick_locus_rep` selects as the representative has 24 aligned
blocks, NONE inside the gene, and one 104,410 bp intron spanning the whole of it. The resulting copy is
154,585 bp for a 23,125 bp gene -- not an over-extended model of NPIPB9 but a different transcript
occupying its coordinates.

The mis-chain filter cannot catch this: it drops introns over 50 kb carried by FEWER THAN 3 reads, and this
one carries 34.

Reported per emitted copy: what fraction of the reads that overlap it contribute no exonic base to it, and
how much of its read support that accounts for. A copy whose support is largely span-only is not a model of
whatever gene lies at its coordinates.
"""
import subprocess
import sys
from collections import defaultdict

BAM, COP, CHROM = sys.argv[1], sys.argv[2], sys.argv[3]
LO = int(sys.argv[4]) if len(sys.argv) > 4 else 0
HI = int(sys.argv[5]) if len(sys.argv) > 5 else 10**9

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 7 or f[3] != CHROM:
        continue
    s, e = int(f[4]), int(f[5])
    if s < LO or e > HI:
        continue
    copies.append({"fam": f[0], "s": s, "e": e, "nx": int(f[6])})
print(f"{len(copies)} emitted copies on {CHROM}:{LO:,}-{HI:,}\n")

rows = []
for c in copies:
    out = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{CHROM}:{c['s']}-{c['e']}"],
                         capture_output=True, text=True).stdout
    tot = span_only = 0
    tot_exonic = 0
    for ln in out.splitlines():
        f = ln.split("\t")
        p, n = int(f[3]) - 1, 0
        inside = 0
        spans = False
        for ch in f[5]:
            if ch.isdigit():
                n = n * 10 + ord(ch) - 48
            else:
                if ch in "M=X":
                    a, b = max(p, c["s"]), min(p + n, c["e"])
                    if b > a:
                        inside += b - a
                    p += n
                elif ch == "N":
                    if p < c["s"] and p + n > c["e"]:
                        spans = True
                    p += n
                elif ch == "D":
                    p += n
                n = 0
        tot += 1
        tot_exonic += inside
        if inside == 0:
            span_only += 1
    if tot:
        rows.append((c, tot, span_only, tot_exonic))

n = len(rows)
frac = sorted(so / t for _, t, so, _ in rows)
print(f"{'reads contributing NO exonic base to the copy they are counted at':<62}")
print(f"  median fraction per copy : {frac[n//2]:.2f}")
print(f"  copies where it exceeds 25% : {sum(1 for x in frac if x>0.25)}/{n} "
      f"({100*sum(1 for x in frac if x>0.25)/n:.0f}%)")
print(f"  copies where it exceeds 50% : {sum(1 for x in frac if x>0.50)}/{n} "
      f"({100*sum(1 for x in frac if x>0.50)/n:.0f}%)")
tot_all = sum(t for _, t, _, _ in rows)
so_all = sum(so for _, _, so, _ in rows)
print(f"  pooled over all copies      : {so_all:,}/{tot_all:,} = {100*so_all/tot_all:.0f}% of read support\n")

big = sorted(rows, key=lambda r: -(r[2] / r[1]))[:10]
print(f"{'copy':<30}{'size':>10}{'reads':>8}{'span-only':>11}{'frac':>7}")
for c, t, so, _ in big:
    print(f"{CHROM + ':' + format(c['s'], ','):<30}{c['e']-c['s']:>10,}{t:>8}{so:>11}{so/t:>7.2f}")

# is span-only support associated with over-sized copies?
sz = [(c["e"] - c["s"], so / t) for c, t, so, _ in rows]
big_c = [f for s, f in sz if s > 50000]
small_c = [f for s, f in sz if s <= 50000]
import statistics as st
if big_c and small_c:
    print(f"\ncopies > 50 kb : median span-only fraction {st.median(big_c):.2f} (n={len(big_c)})")
    print(f"copies <= 50 kb: median span-only fraction {st.median(small_c):.2f} (n={len(small_c)})")
