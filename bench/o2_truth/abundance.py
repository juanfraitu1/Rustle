#!/usr/bin/env python3
"""Per-copy abundance of each tool against the planted reads (register 1970).

    abundance.py PLANT_DIR SCENARIO RUN_DIR      (RUN_DIR = the output folder of bench/o2_families/run_tools.sh)

Reads per planted copy from each tool's own numbers: Rustle `em` of rg.o2.em.tsv (loci of rg.o2.gtf), StringTie `longcov`,
FLAIR fl.isoform.counts.txt (blocks from fl.isoforms.bed), IsoQuant x.discovered_transcript_counts.tsv, bambu `readCount`
(bambu.tsv), and the plain count of primary alignments. A transcript goes to the copy it shares the most exon bases with
(>= 50). Within each family of >= 2 copies: the copy's share of the family's reads against its planted share.
Prints per tool: mean absolute share error (all copies; by unit class: identical to another copy, 1-5 edits, more), and
the share of family reads put on silent copies."""
import sys, os, re, collections

def best_copy(c, e, loci):
    b = (49, None)
    for g, (gc, ge) in loci.items():
        if gc == c and ge[0][0] < e[-1][1] and ge[-1][1] > e[0][0]:
            o = sum(max(0, min(y, q) - max(x, p)) for x, y in e for p, q in ge)
            if o > b[0]: b = (o, g)
    return b[1]

def gtf_tx(path):
    ex = collections.defaultdict(list); meta = {}; at = {}
    for l in open(path):
        f = l.split('\t')
        if len(f) < 9: continue
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if not t: continue
        t = t.group(1)
        if f[2] == 'transcript': at[t] = f[8]
        if f[2] == 'exon': ex[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], re.search(r'gene_id "([^"]+)"', f[8]).group(1))
    return {t: (meta[t][0], sorted(e), meta[t][1], at.get(t, '')) for t, e in ex.items()}

def shares(counts, fam):
    out = {}
    for f, gs in fam.items():
        tot = sum(counts.get(g, 0.0) for g in gs)
        for g in gs: out[g] = counts.get(g, 0.0) / tot if tot > 0 else None
    return out

def tools(d, loci):
    """-> {tool: {copy: reads}}; a tool whose files are missing is left out"""
    T = {}
    def add(name, items):
        c = collections.Counter()
        for chrom, e, v in items:
            g = best_copy(chrom, e, loci)
            if g: c[g] += v
        T[name] = c
    p = d + '/rg.o2.gtf'
    if os.path.exists(p) and os.path.exists(d + '/rg.o2.em.tsv'):
        em = {}
        for l in list(open(d + '/rg.o2.em.tsv'))[1:]: f = l.rstrip('\n').split('\t'); em[f[0]] = float(f[3])
        loc = collections.defaultdict(list); lc = {}
        for t, (c, e, g, a) in gtf_tx(p).items(): loc[g] += e; lc[g] = c
        add('Rustle (EM of the labels)', [(lc[g], sorted(set(e)), em.get(g, 0.0)) for g, e in loc.items()])
    p = d + '/stringtie.gtf'
    if os.path.exists(p): add('StringTie', [(c, e, float(re.search(r'longcov "([^"]+)"', a).group(1)) if 'longcov' in a else 0.0) for c, e, g, a in gtf_tx(p).values()])
    if os.path.exists(d + '/fl.isoform.counts.txt') and os.path.exists(d + '/fl.isoforms.bed'):
        cnt = dict((l.rstrip('\n').rsplit('\t', 1)) for l in open(d + '/fl.isoform.counts.txt')); it = []
        for l in open(d + '/fl.isoforms.bed'):
            f = l.rstrip('\n').split('\t'); s = int(f[1]); sz = [int(x) for x in f[10].strip(',').split(',')]; st = [int(x) for x in f[11].strip(',').split(',')]
            v = cnt.get(f[3]) or next((cnt[k] for k in cnt if k.startswith(f[3]) or f[3].startswith(k)), 0)
            it.append((f[0], [(s + a, s + a + b) for a, b in zip(st, sz)], float(v)))
        add('FLAIR', it)
    p = d + '/iq/x/x.transcript_models.gtf'; q = d + '/iq/x/x.discovered_transcript_counts.tsv'
    if os.path.exists(p) and os.path.exists(q):
        cnt = {l.split('\t')[0]: float(l.split('\t')[1]) for l in list(open(q))[1:] if not l.startswith('__')}
        add('IsoQuant', [(c, e, cnt.get(t, 0.0)) for t, (c, e, g, a) in gtf_tx(p).items()])
    p = d + '/bambu.gtf'; q = d + '/bambu.tsv'
    if os.path.exists(p) and os.path.exists(q):
        R = [l.rstrip('\n').split('\t') for l in open(q)]; i = R[0].index('readCount'); cnt = {r[0]: float(r[i]) if r[i] not in ('NA', '') else 0.0 for r in R[1:]}
        add('bambu', [(c, e, cnt.get(t, 0.0)) for t, (c, e, g, a) in gtf_tx(p).items()])
    if os.path.exists(d + '/x.bam'):
        import pysam
        add('primary alignments, counted', [(r.reference_name, r.get_blocks(), 1.0) for r in pysam.AlignmentFile(d + '/x.bam') if not (r.is_unmapped or r.is_secondary or r.is_supplementary)])
    return T

def main(plant, sc, d):
    plan = [l.rstrip('\n').split('\t') for l in open(plant + '/plan.tsv')]; h = plan[0].index(sc); n = collections.Counter(); famof = {}
    for r in plan[1:]: n[r[1]] += int(r[h]); famof[r[1]] = r[0]
    fam = collections.defaultdict(list)
    for g, f in famof.items(): fam[f].append(g)
    fam = {f: gs for f, gs in fam.items() if len(gs) > 1 and sum(n[g] for g in gs) > 0}
    ex = collections.defaultdict(list); lc = {}
    for l in open(plant + '/loci.gtf'):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); ex[g].append((int(f[3]) - 1, int(f[4]))); lc[g] = f[0]
    loci = {g: (lc[g], sorted(set(e))) for g, e in ex.items()}
    U = [l.rstrip('\n').split('\t') for l in open(plant + '/units.tsv')]; hu = U[0]; size = {r[hu.index('copy')]: int(r[hu.index('unit_size')]) for r in U[1:]}
    near = collections.defaultdict(lambda: 10 ** 9)
    for l in list(open(plant + '/dist.tsv'))[1:]:
        f, a, b, e = l.rstrip('\n').split('\t'); near[a] = min(near[a], int(e)); near[b] = min(near[b], int(e))
    klass = lambda g: 'identical to another copy' if size.get(g, 1) > 1 or near[g] == 0 else '1-5 edits from the nearest' if near[g] <= 5 else 'more than 5 edits'
    truth = shares({g: float(v) for g, v in n.items()}, fam); K = ('identical to another copy', '1-5 edits from the nearest', 'more than 5 edits')
    print(f'## {sc} {d}: {len(fam)} families, {sum(len(v) for v in fam.values())} copies ({sum(1 for gs in fam.values() for g in gs if n[g] == 0)} silent)')
    print('tool\tmean absolute share error\t' + '\t'.join(f'{k} ({sum(1 for gs in fam.values() for g in gs if klass(g) == k)})' for k in K) + '\tshare of family reads on silent copies\tfamilies with no reads from the tool')
    for name, c in tools(d, loci).items():
        s = shares(c, fam); err = collections.defaultdict(list); sil = []; empty = 0
        for f, gs in fam.items():
            if s[gs[0]] is None: empty += 1; continue
            for g in gs: e = abs(s[g] - truth[g]); err['all'].append(e); err[klass(g)].append(e)
            if any(n[g] == 0 for g in gs): sil.append(sum(s[g] for g in gs if n[g] == 0))
        m = lambda v: f'{sum(v) / len(v):.3f}' if v else 'NA'
        print(name, m(err['all']), *[m(err[k]) for k in K], m(sil), empty, sep='\t')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        loci = {'A': ('c', [(0, 100), (200, 300)]), 'B': ('c', [(1000, 1200)])}
        assert best_copy('c', [(50, 100), (200, 260)], loci) == 'A' and best_copy('c', [(1150, 1190)], loci) is None and best_copy('c', [(1000, 1100)], loci) == 'B'
        s = shares({'A': 30.0, 'B': 10.0}, {'f': ['A', 'B'], 'g': ['C', 'D']}); assert s['A'] == 0.75 and s['B'] == 0.25 and s['C'] is None; print('abundance selftest ok')
    else: main(*sys.argv[1:])
