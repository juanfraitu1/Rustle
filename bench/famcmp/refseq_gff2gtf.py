#!/usr/bin/env python3
"""refseq_gff2gtf.py REFSEQ.gff OUT.gtf[.gz]: RefSeq GFF3 -> the truth GTF form the chain scripts read
(one transcript per exon parent; a gene with exons and no transcript child is its own transcript, as in
GGO.GCF_029281585.2_RefSeq.gtf.gz).   --check OLD.gtf.gz NEW.gtf: same set of (contig, strand, exons) models."""
import sys, re, gzip, collections
op = lambda p, m='rt': gzip.open(p, m) if p.endswith('.gz') else open(p, m[0])

def convert(src, dst):
    par = {}; name = {}; ex = collections.defaultdict(list)
    for l in op(src):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if f[2] == 'exon':
            for p in a['Parent'].split(','): ex[(f[0], f[6], p)].append((int(f[3]), int(f[4]), f[1]))
        elif 'ID' in a:
            par[a['ID']] = a.get('Parent', ''); name[a['ID']] = a.get('gene', a.get('Name', ''))
    with op(dst, 'wt') as o:
        for (c, s, t), e in ex.items():
            g = par.get(t, ''); at = f'transcript_id "{t}"; gene_id "{g}"; gene_name "{name.get(g, "")}"'; e.sort()
            o.write(f'{c}\t{e[0][2]}\ttranscript\t{e[0][0]}\t{e[-1][1]}\t.\t{s}\t.\t{at}\n')
            for i, (a_, b, src_) in enumerate(e if s == '+' else e[::-1]):
                o.write(f'{c}\t{src_}\texon\t{a_}\t{b}\t.\t{s}\t.\t{at}; exon_number "{i + 1}";\n')

def models(p):
    m = collections.defaultdict(list)
    for l in op(p):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': m[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    return collections.Counter((k[0], k[1], tuple(sorted(e))) for k, e in m.items())

if __name__ == '__main__':
    if sys.argv[1] == '--check':
        a, b = models(sys.argv[2]), models(sys.argv[3]); print(f'old {sum(a.values())} new {sum(b.values())} only-old {sum((a - b).values())} only-new {sum((b - a).values())}'); sys.exit(a != b)
    convert(sys.argv[1], sys.argv[2])
