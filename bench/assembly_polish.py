#!/usr/bin/env python3
"""Post-assembly precision polish for `copy_assign --assemble-only` GTFs.

Two filters, both using ONLY the emitted `reads "N"` attribute -- no reference, no annotation
(`docs/PREREG_assembly_polish_2026-09-19.md`):

1. Support-aware ISM collapse -- drop a transcript whose intron chain is a contiguous sub-chain of
   another's on the same contig/strand, unless it carries at least `--support-ratio` x the
   container's read support. A mono-exonic transcript inside a multi-exon transcript's span is
   handled the same way.
2. Mono-exonic support floor -- a surviving single-exon transcript has no junction evidence, so it
   must reach the upper quartile of multi-exon read support in the same run. The threshold is
   self-tuning: `--mono-quantile 0.75` of the multi-exon `reads` distribution. `--mono-floor N`
   overrides it with a fixed value; `--mono-quantile 0` disables the filter.

usage: assembly_polish.py IN.gtf OUT.gtf [--support-ratio 1.0] [--mono-quantile 0.75] [--mono-floor N]
"""
import sys, re, collections, argparse

ap = argparse.ArgumentParser()
ap.add_argument("inp"); ap.add_argument("out")
ap.add_argument("--support-ratio", type=float, default=1.0,
                help="keep an ISM fragment when reads(frag) >= ratio * reads(container); 999 = unconditional collapse")
ap.add_argument("--mono-quantile", type=float, default=0.75,
                help="mono-exonic floor = this quantile of multi-exon read support; 0 disables")
ap.add_argument("--mono-floor", type=int, default=None, help="fixed mono-exonic read floor (overrides --mono-quantile)")
ap.add_argument("--no-ism", action="store_true", help="skip filter 1")
a = ap.parse_args()

rows = collections.defaultdict(list); reads = {}
for l in open(a.inp):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    m = re.search(r'transcript_id "([^"]+)"', f[8])
    if not m: continue
    t = m.group(1)
    if f[2] == 'exon':
        rows[t].append((f[0], f[6], int(f[3]) - 1, int(f[4])))
    r = re.search(r'reads "(\d+)"', f[8])
    if r: reads[t] = max(reads.get(t, 0), int(r.group(1)))

chain = {}; span = {}
for t, ex in rows.items():
    ex.sort(key=lambda x: x[2])
    chain[t] = (ex[0][0], ex[0][1], tuple((ex[i][3], ex[i + 1][2]) for i in range(len(ex) - 1)))
    span[t] = (ex[0][2], ex[-1][3])

drop = set()

if not a.no_ism:
    by = collections.defaultdict(list)
    for t, (ch, st, _) in chain.items(): by[(ch, st)].append(t)

    def supported(frag, cont):
        rc = reads.get(cont, 0)
        return rc > 0 and reads.get(frag, 0) >= a.support_ratio * rc

    for key, ts in by.items():
        # deterministic: longest chain first, ties broken by transcript id (the container scan and the
        # mono-exonic host search both depend on this order)
        multi = sorted([t for t in ts if chain[t][2]], key=lambda t: (-len(chain[t][2]), t))
        for x in multi:
            if x in drop: continue
            cx = chain[x][2]
            for y in multi:
                if y == x or y in drop: continue
                cy = chain[y][2]
                if len(cy) >= len(cx): continue
                if any(cx[k:k + len(cy)] == cy for k in range(len(cx) - len(cy) + 1)) and not supported(y, x):
                    drop.add(y)
        for t in sorted(ts):
            if chain[t][2] or t in drop: continue
            s = span[t]
            host = next((m for m in multi if m not in drop and span[m][0] <= s[0] and s[1] <= span[m][1]), None)
            if host is not None and not supported(t, host): drop.add(t)
n_ism = len(drop)

floor = a.mono_floor
if floor is None and a.mono_quantile > 0:
    multi_reads = sorted(reads.get(t, 0) for t in chain if chain[t][2] and t not in drop)
    floor = multi_reads[min(int(a.mono_quantile * len(multi_reads)), len(multi_reads) - 1)] if multi_reads else 0
if floor:
    for t in chain:
        if not chain[t][2] and t not in drop and reads.get(t, 0) < floor: drop.add(t)

with open(a.out, "w") as fo:
    for l in open(a.inp):
        if l.startswith('#'): fo.write(l); continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if m and m.group(1) in drop: continue
        fo.write(l)
print(f"{len(rows)} transcripts -> ISM {n_ism} dropped (ratio {a.support_ratio}) "
      f"-> mono floor {floor} dropped {len(drop) - n_ism} -> {len(rows) - len(drop)} kept", file=sys.stderr)
