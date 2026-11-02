#!/usr/bin/env python3
"""register 1497: drop the transcripts of a GTF whose exons ALL share bases with exons of >= 2 genes of an annotation (same strand),
i.e. the shared constant-exon units of alternative-promoter clusters. usage: drop_shared_units.py IN.gtf ANNOT.gff3 OUT.gtf"""
import collections, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import units_from_annotation as U
gtf, gff, out = sys.argv[1:4]
gix = U.GeneIndex(U.read_genes(gff))
tx = collections.defaultdict(list)
lines = []
for l in open(gtf):
    f = l.split('\t')
    lines.append(l)
    if len(f) > 8 and f[2] == 'exon':
        tx[re.search(r'transcript_id "([^"]+)"', f[8]).group(1)].append((f[0], f[6], int(f[3]), int(f[4])))
drop = {t for t, ex in tx.items() if all(len(gix.hits(c, a - 1, b, s)) >= 2 for c, s, a, b in ex)}
with open(out, 'w') as o:
    for l in lines:
        m = re.search(r'transcript_id "([^"]+)"', l)
        if not (m and m.group(1) in drop):
            o.write(l)
print(f'drop_shared_units: {len(drop)} of {len(tx)} transcripts dropped', file=sys.stderr)
