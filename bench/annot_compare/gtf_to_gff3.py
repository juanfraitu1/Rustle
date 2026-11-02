#!/usr/bin/env python3
"""GTF (exon lines with gene_id / transcript_id / gene_name) -> GFF3 with gene, mRNA and exon records (ID / Parent / Name), the form
bench/units_from_annotation.py reads. usage: gtf_to_gff3.py IN.gtf[.gz] > OUT.gff3"""
import collections, gzip, re, sys
p = sys.argv[1]
fh = gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
G, T, E = {}, {}, []
for l in fh:
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon':
        continue
    g, t = re.search(r'gene_id "([^"]+)"', f[8]), re.search(r'transcript_id "([^"]+)"', f[8])
    if not (g and t):
        continue
    g, t = g.group(1), t.group(1)
    n = re.search(r'gene_name "([^"]+)"', f[8])
    a, b = int(f[3]), int(f[4])
    d = G.setdefault(g, [f[0], a, b, f[6], n.group(1) if n else g])
    d[1], d[2] = min(d[1], a), max(d[2], b)
    x = T.setdefault(t, [g, a, b])
    x[1], x[2] = min(x[1], a), max(x[2], b)
    E.append((t, f[0], a, b, f[6]))
o = sys.stdout
o.write('##gff-version 3\n')
for g, (c, a, b, s, n) in G.items():
    o.write(f'{c}\tconv\tgene\t{a}\t{b}\t.\t{s}\t.\tID={g};Name={n}\n')
for t, (g, a, b) in T.items():
    o.write(f'{G[g][0]}\tconv\tmRNA\t{a}\t{b}\t.\t{G[g][3]}\t.\tID={t};Parent={g}\n')
for k, (t, c, a, b, s) in enumerate(E):
    o.write(f'{c}\tconv\texon\t{a}\t{b}\t.\t{s}\t.\tID=exon{k};Parent={t}\n')
