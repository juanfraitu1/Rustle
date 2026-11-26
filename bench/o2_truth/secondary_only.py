#!/usr/bin/env python3
"""Loci of an assembly that only secondary alignments support, and what speaks for them (register 1971).

    secondary_only.py ASSEMBLY.gtf BAM [--other BAM] [--plant PLANT_DIR SCENARIO] [--truth READS.tsv] [--annot gorilla|human]
                      [--loci FAM.copies.tsv] [--threads N]

Per locus (gene_id, merged exons): primary alignments (any mapping quality) and distinct reads with a secondary alignment
on >= 50 exon bases. primary-supported = >= 1 primary; secondary-only = none and >= 2 such reads.
--other: primaries with mapping quality >= 1 on the locus in another BAM of the same reference. --plant: the planted copy
under the locus is expressed / silent / none. --truth: neutral_truth.py's table, reads whose origin locus overlaps.
--annot: an annotated gene's exons under the locus (env as famtruth). --loci: only the gene_ids of that table (column 13)."""
import sys, os, re, collections, multiprocessing
G = {}

def ov(bl, ex): return sum(max(0, min(b, y) - max(a, x)) for a, b in bl for x, y in ex)

def count(job):
    import pysam
    out = []; B = pysam.AlignmentFile(G['bam']); O = pysam.AlignmentFile(G['other']) if G.get('other') else None
    for g in job:
        c, ex = G['loci'][g]; p = 0; sec = set(); spl = 0
        for r in B.fetch(c, ex[0][0], ex[-1][1]):
            if r.is_supplementary or r.is_unmapped or ov(r.get_blocks(), ex) < 50: continue
            if r.is_secondary:
                if len(sec) < 5000: sec.add(r.query_name)
                spl += 'N' in r.cigarstring
            else: p += 1
        o = None
        if O is not None: o = sum(1 for r in O.fetch(c, ex[0][0], ex[-1][1]) if not (r.is_secondary or r.is_supplementary) and r.mapping_quality >= 1 and ov(r.get_blocks(), ex) >= 50)
        out.append((g, p, len(sec), spl, o))
    return out

def main(argv):
    import pysam
    gtf, bam = argv[:2]; a = argv[2:]; opt = {}
    while a:
        k = a.pop(0)
        if k == '--plant': opt['plant'] = (a.pop(0), a.pop(0))
        else: opt[k[2:]] = a.pop(0)
    ex = collections.defaultdict(list); lc = {}; ntx = collections.Counter(); multi = set(); tex = collections.Counter()
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); ex[g].append((int(f[3]) - 1, int(f[4]))); lc[g] = f[0]; t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); tex[(g, t)] += 1
    for (g, t), n in tex.items():
        if n > 1: multi.add(g)
    def merged(e):
        m = []
        for x, y in sorted(e):
            if m and x <= m[-1][1]: m[-1][1] = max(m[-1][1], y)
            else: m.append([x, y])
        return m
    keep = None
    if 'loci' in opt: keep = {l.split('\t')[12] for l in list(open(opt['loci']))[1:]}
    refs = set(pysam.AlignmentFile(bam).references)
    loci = {g: (lc[g], merged(e)) for g, e in ex.items() if lc[g] in refs and (keep is None or g in keep)}
    G.update(bam=bam, other=opt.get('other'), loci=loci); ids = sorted(loci, key=lambda g: (loci[g][0], loci[g][1][0][0])); thr = int(opt.get('threads', 32)); n = max(1, len(ids) // (thr * 8))
    with multiprocessing.get_context('fork').Pool(thr) as pool: res = pool.map(count, [ids[i:i + n] for i in range(0, len(ids), n)], chunksize=1)
    R = {g: (p, s, spl, o) for r in res for g, p, s, spl, o in r}
    P = [g for g in ids if R[g][0] > 0]; S = [g for g in ids if R[g][0] == 0 and R[g][1] >= 2]
    print(f'loci {len(ids)}: primary-supported {len(P)}, secondary-only (>= 2 reads) {len(S)}, neither {len(ids) - len(P) - len(S)}')
    med = lambda v: sorted(v)[len(v) // 2] if v else 0
    print(f'secondary-only loci: spliced model {sum(g in multi for g in S)} ({sum(g in multi for g in S) / max(1, len(S)):.3f}); reads with a secondary there, median {med([R[g][1] for g in S])}, 2-4 reads {sum(R[g][1] < 5 for g in S)}, 5-19 {sum(5 <= R[g][1] < 20 for g in S)}, 20 or more {sum(R[g][1] >= 20 for g in S)}')
    if opt.get('other'):
        for lab, v in (('primary-supported', P), ('secondary-only', S)): print(f'{lab} loci with a primary alignment (mapping quality >= 1) in the other animal: {sum(R[g][3] > 0 for g in v)} of {len(v)} ({sum(R[g][3] > 0 for g in v) / max(1, len(v)):.3f}); with 5 or more: {sum(R[g][3] >= 5 for g in v)} ({sum(R[g][3] >= 5 for g in v) / max(1, len(v)):.3f})')
    if 'plant' in opt:
        d, sc = opt['plant']; plan = [l.rstrip('\n').split('\t') for l in open(d + '/plan.tsv')]; h = plan[0].index(sc); n = collections.Counter(); pe = collections.defaultdict(list); pc = {}
        for r in plan[1:]: n[r[1]] += int(r[h])
        for l in open(d + '/loci.gtf'):
            f = l.split('\t')
            if len(f) > 8 and f[2] == 'exon': g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); pe[g].append((int(f[3]) - 1, int(f[4]))); pc[g] = f[0]
        def klass(g):
            c, e = loci[g]; b = (49, None)
            for q, qe in pe.items():
                if pc[q] == c:
                    o = ov(e, qe)
                    if o > b[0]: b = (o, q)
            return 'no planted copy' if b[1] is None else 'expressed copy' if n[b[1]] > 0 else 'silent copy'
        for lab, v in (('primary-supported', P), ('secondary-only', S)): k = collections.Counter(klass(g) for g in v); print(f'{lab} loci by planted truth: ' + ', '.join(f'{x} {k[x]}' for x in ('expressed copy', 'silent copy', 'no planted copy')))
    if 'truth' in opt:
        orig = collections.Counter(); cand = collections.Counter(); by = collections.defaultdict(list)
        for g in P + S: by[loci[g][0]].append(g)
        for l in list(open(opt['truth']))[1:]:
            f = l.rstrip('\n').split('\t'); o = set(f[5].split(',')) if f[5] else set()
            for x in f[6].split(';'):
                k, c, se = x.split(':')[:3]; s, e = map(int, se.split('-'))
                for g in by.get(c, ()):
                    if loci[g][1][0][0] < e and loci[g][1][-1][1] > s and ov([(s, e)], loci[g][1]) >= 50: (orig if k in o else cand)[g] += 1
        for lab, v in (('primary-supported', P), ('secondary-only', S)): print(f'{lab} loci in the truth sample: with a read that originates there {sum(orig[g] > 0 for g in v)}, only as another candidate {sum(orig[g] == 0 and cand[g] > 0 for g in v)}, not sampled {sum(orig[g] == 0 and cand[g] == 0 for g in v)}')
        print('   secondary-only loci with an origin read: ' + '; '.join(f'{g} {loci[g][0]}:{loci[g][1][0][0]} origin {orig[g]} candidate {cand[g]} secondary reads {R[g][1]}' for g in S if orig[g] > 0)[:700])
    if 'annot' in opt:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval')); import famtruth as F, bisect
        A = F.load_ggo(F.load_cat()) if opt['annot'] == 'gorilla' else F.load_cat(); ae = collections.defaultdict(list)
        for q, d in A.items():
            for x, y in d['exons']: ae[d['chrom']].append((x, y, d['biotype']))
        for v in ae.values(): v.sort()
        st = {c: [x[0] for x in v] for c, v in ae.items()}; mx = {c: max(y - x for x, y, _ in v) for c, v in ae.items()}
        def bio(g):
            c, e = loci[g]; out = set(); v = ae.get(c, [])
            for x, y in e:
                k = bisect.bisect_left(st.get(c, []), x - mx.get(c, 0))
                while k < len(v) and v[k][0] < y:
                    if v[k][1] > x: out.add(v[k][2])
                    k += 1
            return 'no annotated exon' if not out else 'protein-coding gene' if 'protein_coding' in out else 'pseudogene' if any('pseudo' in b for b in out) else 'other gene'
        for lab, v in (('primary-supported', P), ('secondary-only', S)): k = collections.Counter(bio(g) for g in v); print(f'{lab} loci by annotation: ' + ', '.join(f'{x} {n_} ({n_ / max(1, len(v)):.3f})' for x, n_ in k.most_common()))

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: assert ov([(0, 10), (20, 30)], [(5, 25)]) == 10; print('secondary_only selftest ok')
    else: main(sys.argv[1:])
