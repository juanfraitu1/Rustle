#!/usr/bin/env python3
"""The cluster-level 10-point lead with clusters made from a run's own copy sequences (register 1975).

    cluster_filter_run.py ASSEMBLY.gtf FAM.copies.tsv COPIES_SELF.paf bam   BAM hg TRUTH.tsv GENES.tsv
    cluster_filter_run.py ASSEMBLY.gtf FAM.copies.tsv COPIES_SELF.paf bam   BAM planted PLANT_DIR
    cluster_filter_run.py ASSEMBLY.gtf FAM.copies.tsv COPIES_SELF.paf table NEUTRAL_TRUTH_READS.tsv

COPIES_SELF.paf = PREFIX.fam.copies.fa against itself (minimap2 -cx asm20 -X). Clusters: locus (one locus), 99% and 95%
identity of the copy sequences over >= half of the shorter (single linkage, across family calls); a locus without a
copy sequence is its own cluster. Truth: hg = read -> annotated gene (row 1972) -> the locus under it; planted = the
copy in the read's name; table = the origin loci of neutral_truth.py's table (scores: AS, the table has no ms)."""
import sys, os, re, bisect, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cluster_filter as CF

def load_loci(gtf):
    ex = collections.defaultdict(list); lc = {}
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); ex[g].append((int(f[3]) - 1, int(f[4]))); lc[g] = f[0]
    by = collections.defaultdict(list)
    for g, e in ex.items():
        e = sorted(set(e)); by[lc[g]].append((e[0][0], e[-1][1], e, g))
    for v in by.values(): v.sort(key=lambda x: (x[0], x[1]))
    return by, {c: [x[0] for x in v] for c, v in by.items()}, {c: max(x[1] - x[0] for x in v) for c, v in by.items()}

def locus_of(by, st, mx, c, bl):
    v = by.get(c)
    if not v: return None
    k = bisect.bisect_left(st[c], bl[0][0] - mx[c]); best = (49, None); end = bl[-1][1]
    while k < len(v) and v[k][0] < end:
        s, e, exs, g = v[k]; k += 1
        if e > bl[0][0]:
            o = sum(max(0, min(y, q) - max(x, p)) for x, y in bl for p, q in exs)
            if o > best[0]: best = (o, g)
    return best[1]

def loci_under(by, st, mx, c, bl):
    """every locus sharing >= 50 exon bases with the blocks (two assembled loci can lie on one annotated gene)"""
    v = by.get(c); out = []
    if not v: return out
    k = bisect.bisect_left(st[c], bl[0][0] - mx[c]); end = bl[-1][1]
    while k < len(v) and v[k][0] < end:
        s, e, exs, g = v[k]; k += 1
        if e > bl[0][0] and sum(max(0, min(y, q) - max(x, p)) for x, y in bl for p, q in exs) >= 50: out.append(g)
    return out

def clusters(loci, copies_tsv, paf):
    hdr = {}; rows = [l.rstrip('\n').split('\t') for l in open(copies_tsv)]
    for r in rows[1:]: hdr[f'{r[0]}|{r[1]}'] = r[12]
    out = {'locus': {g: g for g in loci}}
    for thr, lab in ((0.99, '99%'), (0.95, '95%')):
        p, find = CF.uf(); 
        for l in open(paf):
            f = l.split('\t'); a = '|'.join(f[0].split('|')[:2]); b = '|'.join(f[5].split('|')[:2])
            if a != b and a in hdr and b in hdr and int(f[10]) >= 0.5 * min(int(f[1]), int(f[6])) and int(f[9]) / max(1, int(f[10])) >= thr: p[find(hdr[a])] = find(hdr[b])
        out[lab] = {g: find(g) for g in loci}
    return out

def report(title, R, truth, CL, scores):
    size = {k: collections.Counter(v.values()) for k, v in CL.items()}
    print(f'## {title}: ' + '; '.join(f'{k} {len(size[k])} clusters, largest {max(size[k].values())}' for k in CL))
    print('reads\tclusters\tscore\tn\tretained\tshare\twith a truth\ttruth in the retained cluster\tshare\tloci in the retained cluster, median / mean')
    for lab, keep in (('near ties per locus', True), ('all', False)):
        for k, cl in CL.items():
            for si, sn in scores:
                n = ret = wt = ok = 0; sz = []
                for q, al in R.items():
                    top = max(a[si] for a in al)
                    if keep and sum(1 for a in al if a[si] >= top - 10) < 2: continue
                    n += 1; sc = {}
                    for a in al:
                        c = cl[a[0]] if a[0] else None; sc[c] = max(sc.get(c, -10 ** 9), a[si])
                    d = CF.decide(sc)
                    if d is None: continue
                    ret += 1; sz.append(size[k][d]); t = truth.get(q)
                    if t: wt += 1; ok += any(cl[x] == d for x in t)
                sz.sort(); print(f"{lab}\t{k}\t{sn}\t{n}\t{ret}\t{ret / max(1, n):.3f}\t{wt}\t{ok}\t{ok / max(1, wt):.3f}\t{sz[len(sz) // 2] if sz else 0} / {sum(sz) / max(1, len(sz)):.1f}")

def main(gtf, copies, paf, mode, *a):
    by, st, mx = load_loci(gtf); loci = [x[3] for v in by.values() for x in v]; CL = clusters(loci, copies, paf); R = collections.defaultdict(list); truth = {}
    if mode == 'bam':
        import pysam
        bam, tmode = a[0], a[1]
        for r in pysam.AlignmentFile(bam, check_sq=False):
            if r.is_unmapped or r.is_supplementary: continue
            R[r.query_name].append((locus_of(by, st, mx, r.reference_name, r.get_blocks()), r.get_tag('ms'), r.get_tag('AS')))
        if tmode == 'hg':
            gl = {}
            for l in list(open(a[3]))[1:]:
                g, name, fam, c, s, e, stn, ex = l.rstrip('\n').split('\t'); gl[g] = locus_of(by, st, mx, c, [tuple(map(int, x.split('-'))) for x in ex.split(',')])
            for l in open(a[2]):
                q, v, par = l.rstrip('\n').split('\t')
                if v == 'retained' and gl.get(par): truth[q] = [gl[par]]
            print(f'truth paralogs with an assembled locus under them: {sum(1 for v in gl.values() if v)} of {len(gl)}')
        else:
            ex = collections.defaultdict(list); lc = {}
            for l in open(a[2] + '/loci.gtf'):
                f = l.split('\t')
                if len(f) > 8 and f[2] == 'exon': g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); ex[g].append((int(f[3]) - 1, int(f[4]))); lc[g] = f[0]
            cp = {g: locus_of(by, st, mx, lc[g], sorted(set(e))) for g, e in ex.items()}
            for q in R:
                g = q.split('|')[1]
                if cp.get(g): truth[q] = [cp[g]]
        report(f'{tmode} {os.path.basename(os.path.dirname(bam))}', R, truth, CL, ((1, 'ms'), (2, 'AS')))
    else:
        for l in list(open(a[0]))[1:]:
            f = l.rstrip('\n').split('\t'); o = set(f[5].split(',')) if f[5] else set(); t = []
            for x in f[6].split(';'):
                k, c, se, strand, asv = x.split(':')[:5]; s, e = map(int, se.split('-')); g = locus_of(by, st, mx, c, [(s, e)]); R[f[0]].append((g, int(asv[2:]), int(asv[2:])))
                if k in o and g: t.append(g)
            if t: truth[f[0]] = t
        report('neutral-truth table', R, truth, CL, ((2, 'AS'),))

def built(ties, cluster_loci, gtf, mode, *a):
    """score the BUILT step (tools/cluster_verdict.py): TIES.tsv.gz + PREFIX.o2.cluster_loci.tsv against a truth.
    mode hg TRUTH.tsv GENES.tsv | planted PLANT_DIR"""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools')); import cluster_verdict as V
    by, st, mx = load_loci(gtf); cl = {}; size = {}
    for l in list(open(cluster_loci))[1:]: f = l.split('\t'); cl[f[0]] = f[1]; size[f[1]] = int(f[2])
    T = V.read_ties(ties); truth = {}
    if mode == 'hg':
        gl = {}
        for l in list(open(a[1]))[1:]:
            g, name, fam, c, s, e, stn, ex = l.rstrip('\n').split('\t'); gl[g] = loci_under(by, st, mx, c, [tuple(map(int, x.split('-'))) for x in ex.split(',')])
        for l in open(a[0]):
            q, v, par = l.rstrip('\n').split('\t')
            if v == 'retained' and gl.get(par): truth[q] = gl[par]
    else:
        ex = collections.defaultdict(list); lc = {}
        for l in open(a[0] + '/loci.gtf'):
            f = l.split('\t')
            if len(f) > 8 and f[2] == 'exon': g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); ex[g].append((int(f[3]) - 1, int(f[4]))); lc[g] = f[0]
        cp = {g: loci_under(by, st, mx, lc[g], sorted(set(e))) for g, e in ex.items()}
        for q in T:
            g = q.split('|')[1] if '|' in q else None
            if cp.get(g): truth[q] = cp[g]
    n = len(T); dec = wt = ok = 0; sz = []
    for q, rows in T.items():
        v = V.verdict(rows, cl, 10)
        if v is None: continue
        dec += 1; sz.append(size[v]); t = truth.get(q)
        if t: wt += 1; ok += any(cl.get(x) == v for x in t)
    sz.sort(); print(f'built step: tied reads {n}, decided {dec} ({dec / max(1, n):.3f}), with a truth {wt}, right {ok} ({ok / max(1, wt):.3f}), loci in the decided cluster median {sz[len(sz) // 2] if sz else 0} / mean {sum(sz) / max(1, len(sz)):.1f}')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        by = {'c': [(0, 300, [(0, 100), (200, 300)], 'A'), (1000, 1200, [(1000, 1200)], 'B')]}; st = {'c': [0, 1000]}; mx = {'c': 300}
        assert loci_under(by, st, mx, 'c', [(50, 100), (200, 260)]) == ['A'] and loci_under(by, st, mx, 'c', [(0, 1200)]) == ['A', 'B'] and loci_under(by, st, mx, 'c', [(400, 500)]) == []
        assert locus_of(by, st, mx, 'c', [(50, 100), (200, 260)]) == 'A' and locus_of(by, st, mx, 'c', [(1100, 1190)]) == 'B' and locus_of(by, st, mx, 'c', [(400, 500)]) is None and locus_of(by, st, mx, 'd', [(0, 9)]) is None
        print('cluster_filter_run selftest ok')
    elif sys.argv[1] == 'built': built(*sys.argv[2:])
    else: main(*sys.argv[1:])
