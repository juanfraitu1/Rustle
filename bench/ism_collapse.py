#!/usr/bin/env python3
"""Collapse incomplete-splice-match (ISM) fragments out of a transcript GTF.

A transcript is dropped when its intron chain is a CONTIGUOUS SUB-CHAIN of another transcript's chain on the
same contig/strand -- a 5'-truncated piece of something longer that is already emitted (§6p4 localised this
as the 5'UTR truncation signature). Single-exon transcripts are dropped only when contained in a multi-exon
transcript's span.

SUPPORT-AWARE (§6p8): a fragment is KEPT when it carries independent read support, since a well-supported
short chain may be a real shorter isoform rather than a truncation artifact. It is KEPT iff

    reads(fragment) >= --support-abs  (and > 0)          [absolute support]
 or reads(container) > 0 and reads(fragment) >= --support-ratio * reads(container)   [relative support]

and dropped otherwise. The defaults (ratio 999, abs 10**9) make neither clause reachable, i.e. the
unconditional collapse. A transcript with no `reads` attribute counts as 0 and is therefore droppable.

usage: ism_collapse.py IN.gtf OUT.gtf [--support-ratio R] [--support-abs N]
"""
import sys, re, collections, argparse

ap = argparse.ArgumentParser()
ap.add_argument("inp"); ap.add_argument("out")
ap.add_argument("--support-ratio", type=float, default=999.0)
ap.add_argument("--support-abs", type=int, default=10**9)
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

by = collections.defaultdict(list)
for t, (ch, st, c) in chain.items(): by[(ch, st)].append(t)

def keep_for_support(frag, cont):
    rf, rc = reads.get(frag, 0), reads.get(cont, 0)
    if rf > 0 and rf >= a.support_abs: return True
    return rc > 0 and rf >= a.support_ratio * rc

drop = set()
for key, ts in by.items():
    multi = sorted([t for t in ts if chain[t][2]], key=lambda t: -len(chain[t][2]))
    for x in multi:
        if x in drop: continue
        cx = chain[x][2]
        for y in multi:
            if y == x or y in drop: continue
            cy = chain[y][2]
            if len(cy) >= len(cx): continue
            if any(cx[k:k + len(cy)] == cy for k in range(len(cx) - len(cy) + 1)):
                if not keep_for_support(y, x): drop.add(y)
    for t in ts:
        if chain[t][2] or t in drop: continue
        s = span[t]
        host = next((m for m in multi if m not in drop and span[m][0] <= s[0] and s[1] <= span[m][1]), None)
        if host is not None and not keep_for_support(t, host): drop.add(t)

with open(a.out, "w") as fo:
    for l in open(a.inp):
        if l.startswith('#'): fo.write(l); continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if m and m.group(1) in drop: continue
        fo.write(l)
print(f"{len(rows)} transcripts -> dropped {len(drop)} ISM fragments -> {len(rows) - len(drop)} kept "
      f"(ratio {a.support_ratio}, abs {a.support_abs})", file=sys.stderr)
