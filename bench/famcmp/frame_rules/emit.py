"""emit.py TAG DEFAULT_GTF FLOOR1_GTF OUT_GTF : default minus P1 chains plus S1 / S2 chains (multi-exon; single-exon records pass through)."""
import sys, re, collections
n, dp, fp, out = sys.argv[1:5]
ptc = {(f[0], f[1]) for f in (l.rstrip('\n').split('\t') for l in open(n + '.dorf.tsv')) if f[6] == '1'}
add = {tuple(l.rstrip('\n').split('\t')[1:3]) for l in open(n + '.rules.keys')}
def recs(p):
    T = collections.OrderedDict()
    for l in open(p):
        f = l.split('\t')
        if len(f) > 8: T.setdefault(re.search(r'transcript_id "([^"]+)"', f[8]).group(1), []).append(l)
    for t, L in T.items():
        e = sorted((int(x.split('\t')[3]), int(x.split('\t')[4])) for x in L if x.split('\t')[2] == 'exon')
        yield (L[0].split('\t')[0], ','.join('%d-%d' % (e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))), L
done = set()
with open(out, 'w') as w:
    for k, L in recs(dp):
        if k not in ptc: w.writelines(L)
    for k, L in recs(fp):
        if k in add and k not in done: done.add(k); w.writelines(x.replace('transcript_id "', 'transcript_id "S.').replace('gene_id "', 'gene_id "S.') for x in L)
