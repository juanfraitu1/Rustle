#!/usr/bin/env python3
"""Annotated isoforms assembled, by how many reads of their own they have (register 1965).

    thin_isoforms.py ANNOTATION.gff BAM TODAY.gtf ARM=GTF ... [--threads N]

Isoform = (contig, strand, intron chain) of a spliced annotated transcript. Assembled = a transcript of the GTF has exactly
that chain. Own reads = primary alignments with mapping quality >= 1 whose intron chain (the N operations) is exactly the
chain. Prints, per own-read class, the annotated isoforms, those assembled today and by each arm, gained and lost, and the
share of each assembly's spliced transcripts that match an annotated isoform."""
import sys, re, collections, multiprocessing
CLASSES = ('0', '1-2', '3-4', '5-9', '10 or more')
def klass(n): return CLASSES[0] if n == 0 else CLASSES[1] if n < 3 else CLASSES[2] if n < 5 else CLASSES[3] if n < 10 else CLASSES[4]
G = {}

def chains_gtf(path):
    ex = collections.defaultdict(list); meta = {}
    for l in open(path):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); ex[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6])
    out = set()
    for t, e in ex.items():
        e.sort()
        if len(e) > 1: out.add((meta[t][0], meta[t][1], tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))))
    return out

def chains_gff(path):
    ex = collections.defaultdict(list); meta = {}
    for l in open(path):
        if l[0] == '#': continue
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'exon':
            m = re.search(r'Parent=([^;]+)', f[8])
            if m: ex[m.group(1)].append((int(f[3]) - 1, int(f[4]))); meta[m.group(1)] = (f[0], f[6])
    out = set()
    for t, e in ex.items():
        e = sorted(set(e))
        if len(e) > 1: out.add((meta[t][0], meta[t][1], tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))))
    return out

def scan(contig):
    import pysam
    want = G['by'].get(contig)
    if not want: return {}
    n = collections.Counter()
    for r in pysam.AlignmentFile(G['bam']).fetch(contig):
        if r.is_secondary or r.is_supplementary or r.mapping_quality < 1 or 'N' not in r.cigarstring: continue
        c = []; p = r.reference_start
        for o, k in r.cigartuples:
            if o == 3: c.append((p, p + k))
            if o in (0, 2, 3, 7, 8): p += k
        c = tuple(c)
        if c in want: n[c] += 1
    return {(contig, c): v for c, v in n.items()}

def main(gff, bam, today, *arms):
    import pysam
    thr = 16
    if '--threads' in arms: i = arms.index('--threads'); thr = int(arms[i + 1]); arms = arms[:i] + arms[i + 2:]
    ann = chains_gff(gff); by = collections.defaultdict(set)
    for c, s, ch in ann: by[c].add(ch)
    G.update(bam=bam, by=by)
    with multiprocessing.get_context('fork').Pool(thr) as pool: res = pool.map(scan, list(pysam.AlignmentFile(bam).references), chunksize=1)
    own = collections.Counter()
    for d in res: own.update(d)
    kl = {a: klass(own[(a[0], a[2])]) for a in ann}; T = chains_gtf(today); t0 = T & ann
    print(f'annotated isoforms (spliced, distinct chains): {len(ann)}; today: {len(T)} spliced transcripts, {len(t0)} match an annotated isoform ({len(t0) / max(1, len(T)):.4f})')
    print('own reads\tannotated isoforms\tassembled today\t' + '\t'.join(f'{a.split("=")[0]}: assembled (gained / lost)' for a in arms))
    A = {a.split('=', 1)[0]: chains_gtf(a.split('=', 1)[1]) for a in arms}
    for k in CLASSES:
        s = {a for a in ann if kl[a] == k}; row = [k, len(s), len(s & t0)]
        for lab, c in A.items(): m = c & s; row.append(f'{len(m)} (+{len(m - t0)} / -{len((t0 & s) - m)})')
        print(*row, sep='\t')
    thin_missed = {a for a in ann if kl[a] in CLASSES[1:3] and a not in t0}
    for lab, c in A.items():
        m = c & ann; g = m - t0; lost = t0 - m
        print(f'{lab}: {len(c)} spliced transcripts, {len(m)} match ({len(m) / max(1, len(c)):.4f}); gained {len(g)}, lost {len(lost)}; of the gained, thin (1-4 own reads) {sum(kl[a] in CLASSES[1:3] for a in g)}, no own read {sum(kl[a] == CLASSES[0] for a in g)}; thin isoforms missed today {len(thin_missed)}, of them assembled {len(thin_missed & c)} ({len(thin_missed & c) / max(1, len(thin_missed)):.4f})')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        import tempfile, os
        d = tempfile.mkdtemp(); p = os.path.join(d, 'a.gtf'); open(p, 'w').write('c\tx\texon\t1\t10\t.\t+\t.\ttranscript_id "t";\nc\tx\texon\t21\t30\t.\t+\t.\ttranscript_id "t";\nc\tx\texon\t1\t30\t.\t+\t.\ttranscript_id "u";\n')
        q = os.path.join(d, 'a.gff'); open(q, 'w').write('c\tx\texon\t1\t10\t.\t+\t.\tID=e1;Parent=rna1\nc\tx\texon\t21\t30\t.\t+\t.\tID=e2;Parent=rna1\n')
        assert chains_gtf(p) == {('c', '+', ((10, 20),))} == chains_gff(q) and klass(0) == '0' and klass(2) == '1-2' and klass(4) == '3-4' and klass(10) == '10 or more'; print('thin_isoforms selftest ok')
    else: main(*sys.argv[1:])
