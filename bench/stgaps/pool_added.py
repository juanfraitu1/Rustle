#!/usr/bin/env python3
"""Pool-added transcripts: multi-exon transcripts of the pooled run whose (contig, strand, intron chain) is absent from
the primary-only run. Annotation-free. usage: pool_added.py POOLED.gtf PRIMARY.gtf OUT.tsv (tid contig strand chain)"""
import re, sys
from collections import defaultdict
RT = re.compile(r'transcript_id "([^"]+)"')
def chains(p):
    ex, key = defaultdict(list), {}
    for l in open(p):
        f = l.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        t = RT.search(f[8]).group(1)
        ex[t].append((int(f[3]) - 1, int(f[4])))
        key[t] = (f[0], f[6])
    out = {}
    for t, e in ex.items():
        if len(e) < 2:
            continue
        e.sort()
        out[t] = (key[t][0], key[t][1], tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)))
    return out
cur, prim = chains(sys.argv[1]), chains(sys.argv[2])
have = set(prim.values())
n = 0
with open(sys.argv[3], "w") as fo:
    for t, k in sorted(cur.items()):
        if k not in have:
            n += 1
            fo.write(f"{t}\t{k[0]}\t{k[1]}\t{','.join(f'{a}-{b}' for a, b in k[2])}\n")
print(f"pooled multi-exon {len(cur)}, primary-only {len(prim)}, pool-added {n}")
