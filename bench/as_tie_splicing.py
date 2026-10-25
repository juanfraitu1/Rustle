#!/usr/bin/env python3
"""as_tie_splicing.py — is the AS gap between a read's primary and its secondaries a copy signal? (2026-10-24, LEDGER_O2, register 1403)

Input: per region group g, two line-aligned files from the same BAM pass (no unmapped, no supplementary records):
  samtools view -@ 4 -F 0x804 --keep-tag AS,NM,s1,ms BAM <contigs> | cut -f1-5,12- > g.tsv         (name flag chrom pos MAPQ tags)
  samtools view -@ 4 -F 0x804 BAM <contigs> | cut -f1,2,6 | awk '{n=gsub(/N/,"",$3); print $1"\t"$2"\t"n}' > g.intr   (introns)
usage: as_tie_splicing.py g1 [g2 ...]   (prefixes; prints part 1, the AS-gap classes, and part 2, the splice-aware windows)
"""
import collections
import sys

PEN = 32   # splice:hq long-gap (intron) open cost (-O 2,32 -E 1,0): the splice-neutral AS adds it back per intron

prim, secs = {}, collections.defaultdict(list)
for g in sys.argv[1:]:
    for ln, li in zip(open(f"{g}.tsv"), open(f"{g}.intr")):
        f = ln.rstrip("\n").split("\t")
        n_, fl, ni = li.rstrip("\n").split("\t")
        assert n_ == f[0] and fl == f[1], "the .tsv and .intr files are not line-aligned"
        d = {t[:2]: int(t[5:]) for t in f[5:]}
        rec = (d["AS"], d["NM"], d.get("s1", -1), int(f[4]), int(ni))     # AS NM s1 MAPQ introns
        if int(f[1]) & 256:
            secs[f[0]].append(rec)
        else:
            prim[f[0]] = rec

# part 1: the best secondary (by AS, then fewest NM) against the primary
T = collections.defaultdict(collections.Counter)
c = collections.Counter()
for n, ss in secs.items():
    p = prim.get(n)
    if p is None:
        continue
    c["reads with a secondary"] += 1
    b = max(ss, key=lambda r: (r[0], -r[1]))
    d = b[0] - p[0]
    if d >= 0:
        k = "exact tie (0)" if d == 0 else ("higher by 1-9" if d < 10 else "higher by >= 10")
        t = T[k]
        t["reads"] += 1
        t["secondary has fewer introns"] += b[4] < p[4]
        t["secondary unspliced, primary spliced"] += b[4] == 0 and p[4] > 0
        t["secondary has more edits (NM)"] += b[1] > p[1]
        t["same introns and fewer NM"] += b[4] == p[4] and b[1] < p[1]
        t["primary MAPQ >= 20"] += p[3] >= 20
    # part 2: does SOME secondary fall in each window
    windows = (("A raw AS within 10 (or higher)", [s for s in ss if s[0] >= p[0] - 10]),
               ("B splice-neutral AS within 10", [s for s in ss if s[0] + PEN * s[4] >= p[0] + PEN * p[4] - 10]),
               ("C same intron count + AS within 10", [s for s in ss if s[4] == p[4] and abs(s[0] - p[0]) <= 10]))
    for k, v in windows:
        if v:
            c[k] += 1
            c[k + " | a candidate with > primary NM + 5"] += any(s[1] > p[1] + 5 for s in v)
            c[k + " | primary MAPQ >= 20"] += p[3] >= 20
print(f"primaries {len(prim):,}")
for k in ("exact tie (0)", "higher by 1-9", "higher by >= 10"):
    t = T[k]
    print(f"== best secondary AS - primary AS {k}: {t['reads']:,} reads")
    for kk, v in t.items():
        if kk != "reads":
            print(f"   {kk}: {v:,} ({v / max(1, t['reads']):.1%})")
for k, v in c.items():
    print(f"{k}: {v:,}")
