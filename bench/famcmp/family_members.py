#!/usr/bin/env python3
"""Are the annotated copies of one gene family, and their isoforms, assembled? Descriptive, per gene.

  family_members.py ANNOT.gff3[.gz] NAME_REGEX BAM COPIES.tsv NAME=ASSEMBLY.gtf ...
  family_members.py --isoforms ANNOT.gff3[.gz] NAME_REGEX BAM COPIES.tsv ASSEMBLY.gtf OUT_PREFIX [THREADS]
  family_members.py --selftest
Per annotated gene whose gene_name matches NAME_REGEX: annotated isoforms (distinct intron chains), reads whose PRIMARY
record overlaps its exons / reads with ANY record there (secondaries included), and per assembly: models overlapping its
exons on its strand, annotated intron chains matched exactly, and the family (COPIES.tsv) holding a copy on the gene.

--isoforms (2026-11-13, descriptive): every annotated spliced isoform of every gene matching NAME_REGEX that carries a copy
of COPIES.tsv on its strand, in the first class that fits (reads must carry the chain exactly, no junction tolerance):
  exact | near10 (a model within 10 bp at every junction) | longer_fused / longer_own (a model contains the chain and
  does / does not touch exons of another annotated gene on that strand) | noread_part / noread (no record at the locus
  carries the full chain; a model is / is not a contiguous part of it) | secondary_only (full-chain records, none primary)
  | below_floor (1-2 full-chain primaries) | partial_model / missed (>= 3 full-chain primaries; a model is a part / none).
OUT_PREFIX.isoforms.tsv (one row per isoform) and OUT_PREFIX.summary.txt (classes overall and per family)."""
import sys, re, gzip, collections

def chain(exons):
    e = sorted(exons); return tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))

def ov(a, b):
    return any(s < y and x < e for s, e in a for x, y in b)

def read_gtf(path):
    tx = collections.defaultdict(list); meta = {}
    for l in open(path):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            tx[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6])
    return [(meta[t][0], meta[t][1], sorted(e)) for t, e in tx.items()]

def main(annot, rx, bam, copies, arms):
    import pysam
    genes, tgene, tex = {}, {}, collections.defaultdict(list)
    for l in (gzip.open(annot, 'rt') if annot.endswith('.gz') else open(annot)):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if f[2] == 'gene' and re.search(rx, a.get('gene_name', '')): genes[a['ID']] = (a['gene_name'], f[0], f[6], a.get('gene_biotype', ''))
        elif f[2] == 'transcript': tgene[a['ID']] = a['Parent']
        elif f[2] == 'exon' and tgene.get(a['Parent']) in genes: tex[a['Parent']].append((int(f[3]) - 1, int(f[4])))
    iso = collections.defaultdict(set); gex = collections.defaultdict(set)
    for t, e in tex.items(): iso[tgene[t]].add(chain(e)); gex[tgene[t]].update(e)
    fam = [(r[0], r[3], int(r[4]), int(r[5]), r[7]) for r in (l.split('\t') for l in list(open(copies))[1:])]
    models = {n: read_gtf(p) for n, p in arms}
    B = pysam.AlignmentFile(bam)
    print('gene\tbiotype\tlocus\tisoforms\tprimary_reads\tany_reads\tfamily\t' + '\t'.join(f'{n}_models\t{n}_isoforms_exact' for n, _ in arms))
    tot = collections.Counter()
    for g, (name, c, st, bt) in sorted(genes.items(), key=lambda kv: kv[1]):
        ex = sorted(gex[g]); lo, hi = ex[0][0], ex[-1][1]; prim, anyr = set(), set()
        for r in B.fetch(c, lo, hi):
            if r.is_supplementary or not ov(r.get_blocks(), ex): continue
            anyr.add(r.query_name)
            if not r.is_secondary: prim.add(r.query_name)
        fams = sorted({f for f, fc, s, e, fs in fam if fc == c and fs == st and ov([(s, e)], ex)})
        row = [name, bt, f'{c}:{lo}-{hi}{st}', len(iso[g]), len(prim), len(anyr), ','.join(fams) or '-']
        tot['genes'] += 1; tot['isoforms'] += len(iso[g]); tot['in_family'] += bool(fams); tot['expressed'] += len(prim) >= 2; tot['any2'] += len(anyr) >= 2
        for n, _ in arms:
            hit = [m for m in models[n] if m[0] == c and m[1] == st and ov(m[2], ex)]
            exact = len({chain(m[2]) for m in hit} & {i for i in iso[g] if i})
            row += [len(hit), exact]; tot[n + '_genes'] += bool(hit); tot[n + '_iso'] += exact
        print('\t'.join(map(str, row)))
    print('# ' + ', '.join(f'{k} {v}' for k, v in tot.items()))

def load_annot(annot, rx):
    """genes {id: (name, chrom, strand)}, isoforms {gene: {chain}}, exons {gene: [exon]} from a CAT, Ensembl or RefSeq GFF3."""
    genes, par, tex = {}, {}, collections.defaultdict(list)
    for l in (gzip.open(annot, 'rt') if annot.endswith('.gz') else open(annot)):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if f[2] in ('gene', 'pseudogene'): genes[a['ID']] = (a.get('gene_name') or a.get('Name') or a['ID'], f[0], f[6])
        elif f[2] == 'exon': tex[a['Parent']].append((int(f[3]) - 1, int(f[4])))
        elif 'Parent' in a and 'ID' in a: par[a['ID']] = a['Parent']
    iso = collections.defaultdict(set); gex = collections.defaultdict(set)
    for t, e in tex.items():
        g = par.get(t, t)   # RefSeq pseudogene exons hang off the gene itself
        if g in genes and re.search(rx, genes[g][0]): iso[g].add(chain(e)); gex[g].update(e)
    return {g: genes[g] for g in gex}, iso, {g: sorted(e) for g, e in gex.items()}

def introns(r):
    p = r.reference_start; out = []
    for op, n in r.cigartuples:
        if op == 3: out.append((p, p + n))
        if op in (0, 2, 3, 7, 8): p += n
    return tuple(out)

sub = lambda a, b: len(a) >= 1 and any(b[i:i + len(a)] == a for i in range(len(b) - len(a) + 1))
near = lambda a, b: len(a) == len(b) and all(abs(x[0] - y[0]) <= 10 and abs(x[1] - y[1]) <= 10 for x, y in zip(a, b))

def classify(i, mch, fused, rp, ra):
    """mch: chains of the models on the gene; fused: those of models that also touch another gene; rp / ra: read chain counts."""
    if i in mch: return 'exact'
    if any(near(i, m) for m in mch): return 'near10'
    if any(sub(i, m) for m in mch): return 'longer_own' if any(sub(i, m) for m in mch if m not in fused) else 'longer_fused'
    part = any(sub(m, i) for m in mch)
    if ra[i] == 0: return 'noread_part' if part else 'noread'
    if rp[i] == 0: return 'secondary_only'
    if rp[i] < 3: return 'below_floor'
    return 'partial_model' if part else 'missed'

def _chrom(job):
    import pysam
    bam, c, items, models, others = job; B = pysam.AlignmentFile(bam); out = []
    for g, name, st, ex, isos, fams in items:
        lo, hi = ex[0][0], ex[-1][1]; rp, ra = collections.Counter(), collections.Counter()
        want = {i for i in isos if i}
        if want and c in B.references:
            for r in B.fetch(c, lo, hi):
                if r.is_supplementary or not r.cigartuples: continue
                ch = introns(r)
                if ch in want:
                    ra[ch] += 1
                    if not r.is_secondary: rp[ch] += 1
        hit = [m for m in models if m[0] == st and m[1][0][0] < hi and lo < m[1][-1][1] and ov(m[1], ex)]
        mch = {chain(m[1]) for m in hit}
        oth = [e for h, s2, e in others if h != g and s2 == st and e[0][0] < hi + 2000000 and lo - 2000000 < e[-1][1]]
        fused = {chain(m[1]) for m in hit if any(ov(m[1], e) for e in oth)}
        for i in sorted(want): out.append((name, f'{c}:{lo}-{hi}{st}', fams, len(i), rp[i], ra[i], classify(i, mch, fused, rp, ra)))
    return out

def isoforms(annot, rx, bam, copies, gtf, out, threads=8):
    import multiprocessing
    genes, iso, gex = load_annot(annot, rx)
    fam = collections.defaultdict(list)
    for r in (l.split('\t') for l in list(open(copies))[1:]): fam[(r[3], r[7])].append((int(r[4]), int(r[5]), r[0]))
    models = collections.defaultdict(list)
    for c, st, e in read_gtf(gtf): models[c].append((st, e))
    bych = collections.defaultdict(list); others = collections.defaultdict(list)
    for g, (name, c, st) in genes.items():
        others[c].append((g, st, gex[g]))
        fs = sorted({f for s, e, f in fam.get((c, st), ()) if ov([(s, e)], gex[g])})
        if fs: bych[c].append((g, name, st, gex[g], iso[g], ','.join(fs)))
    jobs = [(bam, c, v, models[c], others[c]) for c, v in sorted(bych.items(), key=lambda kv: -len(kv[1]))]
    with multiprocessing.Pool(int(threads)) as P: rows = [r for part in P.imap_unordered(_chrom, jobs) for r in part]
    rows.sort()
    with open(out + '.isoforms.tsv', 'w') as w:
        w.write('gene\tlocus\tfamilies\tintrons\tfull_primary\tfull_any\tclass\n')
        for r in rows: w.write('\t'.join(map(str, r)) + '\n')
    tot = collections.Counter(r[6] for r in rows); n = sum(tot.values())
    perfam = collections.defaultdict(collections.Counter)
    for r in rows: perfam[r[2].split(',')[0]][r[6]] += 1
    with open(out + '.summary.txt', 'w') as w:
        w.write(f'genes in a family {sum(len(v) for v in bych.values())}, spliced isoforms {n}\n')
        for k, v in tot.most_common(): w.write(f'{k}\t{v}\t{v / max(1, n):.3f}\n')
        w.write('-- families with the most secondary_only + missed + partial_model + longer_fused\n')
        key = lambda c: c['secondary_only'] + c['missed'] + c['partial_model'] + c['longer_fused']
        for f, c in sorted(perfam.items(), key=lambda kv: -key(kv[1]))[:25]: w.write(f + '\t' + ' '.join(f'{k}={v}' for k, v in c.most_common()) + '\n')
    print(open(out + '.summary.txt').read())

def selftest():
    C = collections.Counter; a = ((5, 10), (20, 30))
    assert classify(a, {a}, set(), C(), C()) == 'exact'
    assert classify(a, {((6, 10), (20, 31))}, set(), C(), C()) == 'near10'
    big = ((1, 2),) + a
    assert classify(a, {big}, {big}, C(), C()) == 'longer_fused' and classify(a, {big}, set(), C(), C()) == 'longer_own'
    assert classify(a, {a[:1]}, set(), C(), C()) == 'noread_part' and classify(a, set(), set(), C(), C()) == 'noread'
    assert classify(a, set(), set(), C(), C({a: 4})) == 'secondary_only'
    assert classify(a, set(), set(), C({a: 2}), C({a: 2})) == 'below_floor'
    assert classify(a, set(), set(), C({a: 3}), C({a: 3})) == 'missed' and classify(a, {a[1:]}, set(), C({a: 3}), C({a: 3})) == 'partial_model'

    assert chain([(10, 20), (0, 5)]) == ((5, 10),) and chain([(0, 5)]) == ()
    assert ov([(0, 5)], [(4, 9)]) and not ov([(0, 5)], [(5, 9)])
    print('selftest ok')

if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else isoforms(*sys.argv[2:]) if sys.argv[1] == '--isoforms' else main(*sys.argv[1:5], [a.split('=', 1) for a in sys.argv[5:]])
