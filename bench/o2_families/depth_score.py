#!/usr/bin/env python3
"""Score downsampled assemblies against the full-depth primaries-only assembly (register 1964).

    depth_score.py TRUTH.gtf REFSEQ.gff primaries=A.gtf pool=B.gtf [mode=C.gtf ...]

A spliced intron chain (contig, strand, introns) is confirmed when TRUTH.gtf has it. The first arm is the reference:
for every other arm, added = its chains the reference lacks, lost = the reverse. Prints one row per arm."""
import sys, re, collections
def chains(path, gff=False):
    ex = collections.defaultdict(list)
    for l in open(path):
        if l[0] == '#': continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        t = re.search(r'Parent=([^;]+)' if gff else r'transcript_id "([^"]+)"', f[8])
        if t: ex[(t.group(1), f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
    out = set()
    for (t, c, s), e in ex.items():
        e.sort()
        if len(e) > 1: out.add((c, s, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))))
    return out
truth = chains(sys.argv[1]); ref = chains(sys.argv[2], True); arms = [a.split('=', 1) for a in sys.argv[3:]]
assert truth and ref, 'empty truth or annotation'
base = chains(arms[0][1]); bc = len(base & truth)
print('arm\tchains\tconfirmed\tprecision\tadded\tadded_confirmed\tadded_confirmed_share\tadded_refseq_exact\tlost\tlost_confirmed\tgain_share_of_reference_confirmed')
for name, path in arms:
    c = chains(path); add = c - base; lost = base - c; ac = len(add & truth)
    print(name, len(c), len(c & truth), f'{len(c & truth) / max(1, len(c)):.3f}', len(add), ac, f'{ac / max(1, len(add)):.3f}', len(add & ref), len(lost), len(lost & truth), f'{ac / max(1, bc):.4f}', sep='\t')
