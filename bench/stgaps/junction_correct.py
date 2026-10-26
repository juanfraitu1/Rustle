#!/usr/bin/env python3
"""Read-level junction correction for PREREG_junction_correct_2026-10-26.md: StringTie's long-read rule (ST) and its
sequence-checked form (SEQ). Writes a corrected BAM of one contig; the assembler then runs unchanged on it.

usage: junction_correct.py st|seq|none BAM FASTA CONTIG OUT.bam

Sites = (position, strand) of a junction's LEFT end (intron start, 0-based) and RIGHT end (intron end, exclusive),
support = primary spliced reads (flag & 2308 == 0) with an intron there. Strand = the junction's motif (GT-AG, GC-AG,
AT-AC = +; CT-AC, CT-GC, GT-AT = -), else the read's minimap2 `ts` tag (flipped for a reverse alignment), else 0.
ST: each site moves to the same-strand site within < 25 bp with the largest support, if support_t * 0.9 > support_own
(ties: nearer, then left); one step. SEQ: ST's move kept only if the read's bases around the junction (<= 40 bp each
side) have edit distance (edlib NW) to the new spliced genome <= that to the old.
A move is encoded in the CIGAR (left +d: D d before N; left -k: I k; right +e: I e after N; right -k: D k), so exon
coordinates change and query sequence / alignment start do not. Skipped (counted) when the ops flanking the N are not
M/=/X long enough for the shift.
"""
import bisect
import sys
from collections import Counter, defaultdict

import pysam

W, TOL, F = 25, 0.9, 40
SKIP = 0x4 | 0x100 | 0x800
MTYPE = {0, 7, 8}  # M = X
QCONS = {0, 1, 4, 7, 8}
RCONS = {0, 2, 3, 7, 8}
PLUS = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
MINUS = {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}

mode, bam_path, fa_path, contig, out_path = sys.argv[1:6]
if mode not in ("st", "seq", "none"):
    sys.exit("mode must be st, seq or none")
if mode == "seq":
    import edlib
g = pysam.FastaFile(fa_path).fetch(contig).upper()


def strand_of(L, R, rec):
    m = (g[L:L + 2], g[R - 2:R])
    if m in PLUS:
        return 1
    if m in MINUS:
        return -1
    if rec.has_tag("ts"):
        s = 1 if rec.get_tag("ts") == "+" else -1
        return -s if rec.is_reverse else s
    return 0


def junctions(rec):
    """[(op index of N, L, R, query offset at N)]."""
    out, r, q = [], rec.reference_start, 0
    for i, (op, n) in enumerate(rec.cigartuples):
        if op == 3:
            out.append((i, r, r + n, q))
        if op in RCONS:
            r += n
        if op in QCONS:
            q += n
    return out


def eligible(rec):
    return not (rec.flag & SKIP) and rec.cigartuples is not None and any(op == 3 for op, _ in rec.cigartuples)


bam = pysam.AlignmentFile(bam_path, "rb")
left, right = Counter(), Counter()
for rec in bam.fetch(contig):
    if rec.flag & 0x200 or not eligible(rec):
        continue
    for _, L, R, _ in junctions(rec):
        s = strand_of(L, R, rec)
        left[(L, s)] += 1
        right[(R, s)] += 1


def targets(sup):
    by = defaultdict(list)
    for (p, s) in sup:
        by[s].append(p)
    tgt = {}
    for s, ps in by.items():
        ps.sort()
        for p in ps:
            own = sup[(p, s)]
            best = None
            for q in ps[bisect.bisect_right(ps, p - W):bisect.bisect_left(ps, p + W)]:
                if q == p or sup[(q, s)] * TOL <= own:
                    continue
                k = (-sup[(q, s)], abs(q - p), q)
                if best is None or k < best[0]:
                    best = (k, q)
            if best is not None:
                tgt[(p, s)] = best[1]
    return tgt


tl, tr = targets(left), targets(right)
cnt = Counter()


def edit(a, b):
    return edlib.align(a, b, mode="NW")["editDistance"]


out = pysam.AlignmentFile(out_path, "wb", template=bam)
for rec in bam.fetch(contig):
    if mode == "none" or rec.flag & 0x200 or not eligible(rec):
        out.write(rec)
        continue
    ops = [list(x) for x in rec.cigartuples]
    edits = {}  # N op index -> (d, e)
    for i, L, R, q in junctions(rec):
        s = strand_of(L, R, rec)
        cnt["junctions"] += 1
        L2, R2 = tl.get((L, s), L), tr.get((R, s), R)
        d, e = L2 - L, R2 - R
        if d == 0 and e == 0:
            continue
        cnt["candidates"] += 1
        if i == 0 or i + 1 >= len(ops) or ops[i - 1][0] not in MTYPE or ops[i + 1][0] not in MTYPE \
                or (d < 0 and ops[i - 1][1] <= -d) or (e > 0 and ops[i + 1][1] <= e) or R2 - L2 <= 0:
            cnt["skipped_flank"] += 1
            continue
        if mode == "seq":
            fl = min(F, ops[i - 1][1], q)
            fr = min(F, ops[i + 1][1], len(rec.query_sequence) - q)
            if fl <= max(0, -d) or fr <= max(0, e) or fl <= abs(d) or fr <= abs(e):
                cnt["skipped_flank"] += 1
                continue
            rd = rec.query_sequence[q - fl:q + fr].upper()
            old = g[L - fl:L] + g[R:R + fr]
            new = g[L - fl:L2] + g[R2:R + fr]
            if edit(rd, new) > edit(rd, old):
                cnt["rejected_seq"] += 1
                continue
        edits[i] = (d, e)
        cnt["moved"] += 1
    if edits:
        new = []
        for i, (op, n) in enumerate(ops):
            if i in edits:
                d, e = edits[i]
                if d > 0:
                    new.append((2, d))
                elif d < 0:
                    new[-1] = (new[-1][0], new[-1][1] + d)
                    new.append((1, -d))
                new.append((3, n + e - d))
                if e > 0:
                    new.append((1, e))
                elif e < 0:
                    new.append((2, -e))
                if e > 0:
                    ops[i + 1][1] -= e
            else:
                new.append((op, n))
        if any(n <= 0 for _, n in new):  # two moves ate one short exon: leave the read as aligned
            cnt["reverted_short_exon"] += 1
            cnt["moved"] -= len(edits)
            out.write(rec)
            continue
        rec.cigartuples = new
        cnt["reads_changed"] += 1
    out.write(rec)
out.close()
pysam.index(out_path)
print(f"[junction_correct] {mode} {contig}: left sites {len(left)} ({len(tl)} with a target), right sites {len(right)} "
      f"({len(tr)} with a target); {dict(cnt)}", file=sys.stderr)
