#!/usr/bin/env python3
"""O2 family report of one run (driver stage `o2`, the O2 stage of `all` since 2026-11-22, register 1903).

    o2_family_report.py PREFIX.gtf PREFIX.fam.copies.tsv OUT_PREFIX [--eq PREFIX.eq.loci.tsv]

What the assembly reports for each gene family of `families`: its copies (loci), the copies with a spliced transcript,
transcripts, distinct intron chains (isoforms) and reads. With --eq (stage `eqclasses`, or `multicopy`'s
PREFIX.multicopy.eq.loci.tsv) also the reads placed on one copy alone, the tie groups and the range of copies the reads
can show expressed (copies with >= 2 placed reads .. copies with any read). No read is assigned to a copy.
Writes OUT_PREFIX.families.tsv (one row per family) and OUT_PREFIX.copies.tsv (one row per copy).
--sequences PREFIX.o2.sequences.tsv (register 1926, bench/o2_families/recurrence_test.py): adds the column tied_sequences =
the distinct sequences (per isoform) the reads show among copies they are exactly tied between, counted at the copy the
tied reads were read at; in the family row their sum. A count of sequences, not of copies: two alleles of one copy are two.
The same question for several tools, with planted truth: bench/o2_families/ (docs/O2_FAMILY_BENCHMARK.md).
A copy named GENE.pK (an extra representative of a locus, added by the families stage's `famreps` step) is a group of
the transcripts of GENE: each transcript of GENE is counted at the piece whose representative exons it shares the most
bases with (as tools/locus_reps.py groups them). Tie groups and tied sequences are per locus, so they are reported at
the first piece (.p1) only. Before 2026-11-25 such copies were reported with 0 transcripts (54 of 1,345 on OR6737).
--selftest runs the pure parts."""
import sys, os, re, csv, collections

def load_gtf(path):
    """-> {gene_id: {'tx': n, 'spliced': n, 'chains': set, 'reads': float, 'name': str}}"""
    ex = collections.defaultdict(list); meta = {}
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]); g = re.search(r'gene_id "([^"]+)"', f[8])
        if not t or not g: continue
        if f[2] == 'exon': ex[t.group(1)].append((int(f[3]), int(f[4]))); meta.setdefault(t.group(1), [g.group(1), 0.0, '.'])
        else:
            r = re.search(r' reads "([^"]+)"', ' ' + f[8]); n = re.search(r'ref_gene_name "([^"]+)"', f[8])
            m = meta.setdefault(t.group(1), [g.group(1), 0.0, '.'])
            if r: m[1] = float(r.group(1))
            if n: m[2] = n.group(1)
    out = collections.defaultdict(new_gene); TX.clear()
    for t, e in ex.items():
        g, r, n = meta[t]; e.sort(); add_tx(out[g], e, r, n); TX[g].append((e, r, n))
    return out

TX = collections.defaultdict(list)      # gene -> [(exons, reads, name)] of the GTF last loaded, for the pieces of a locus
def new_gene(): return {'tx': 0, 'spliced': 0, 'chains': set(), 'reads': 0.0, 'name': '.'}
def add_tx(d, e, r, n):
    d['tx'] += 1; d['reads'] += r
    if n != '.': d['name'] = n
    if len(e) > 1: d['spliced'] += 1; d['chains'].add(tuple((a[1], b[0]) for a, b in zip(e, e[1:])))

PIECE = re.compile(r'\.p(\d+)$')
def pieces(genes, copies):
    """genes gets an entry for every copy GENE.pK: the transcripts of GENE that share the most exon bases with that piece"""
    by = collections.defaultdict(dict)
    for r in copies:
        if PIECE.search(r['gene_id']) and r['gene_id'] not in genes and r.get('exons'):
            by[PIECE.sub('', r['gene_id'])][r['gene_id']] = [(int(a) + 1, int(b)) for a, b in (x.split('-') for x in r['exons'].split(','))]   # the table is 0-based half-open
    shared = lambda a, b: sum(max(0, min(e, e2) - max(s, s2) + 1) for s, e in a for s2, e2 in b)
    for parent, ps in by.items():
        for e, r, n in TX.get(parent, ()):
            sh = {k: shared(e, x) for k, x in ps.items()}; k = max(sorted(sh), key=lambda k: sh[k])
            if sh[k] > 0: add_tx(genes[k], e, r, n)

def report(genes, copies, eq=None, sq=None):
    """copies: rows of fam.copies.tsv (dicts). -> family rows, copy rows"""
    fam = collections.defaultdict(list); seen = set()
    for r in copies:
        k = (r['family_id'], r['gene_id'])
        if k not in seen: seen.add(k); fam[r['family_id']].append(r)
    F, C = [], []
    for f in sorted(fam, key=lambda f: (-len(fam[f]), f)):
        nsq = 0; tx = sp = iso = 0; reads = 0.0; chroms = set(); names = set(); placed = mx = 0; groups = set(); lo = hi = 0
        for r in fam[f]:
            g = genes.get(r['gene_id'], {'tx': 0, 'spliced': 0, 'chains': set(), 'reads': 0.0, 'name': '.'})
            first = PIECE.sub('', r['gene_id']) if r['gene_id'].endswith('.p1') else None      # per-locus tables: at the first piece only
            e = (eq or {}).get(r['gene_id']) or (eq or {}).get(first)
            tx += g['tx']; sp += g['spliced'] > 0; iso += len(g['chains']); reads += g['reads']; chroms.add(r['chrom'])
            if g['name'] != '.': names.add(g['name'])
            row = [f, r['gene_id'], r['chrom'], r.get('locus_start', r['start']), r.get('locus_end', r['end']), r['strand'], g['name'], g['tx'], len(g['chains']), f"{g['reads']:g}"]
            if eq is not None:
                p, m = (int(e['placed']), int(e['max'])) if e else (0, 0)
                placed += p; mx += m; lo += p >= 2; hi += m > 0
                if e and e['group'] != '.': groups.add(e['group'])
                row += [p, m, e['em'] if e else '0.0', e['group'] if e else '.']
            if sq is not None: k = sq.get(r['gene_id'], 0) or sq.get(first, 0); row.append(k); nsq += k
            C.append(row)
        row = [f, len(fam[f]), sp, tx, iso, f'{reads:g}', len(chroms), ','.join(sorted(names)) or '.']
        if eq is not None: row += [placed, len(groups), max(1, lo) if hi else 0, hi]
        if sq is not None: row.append(nsq)
        F.append(row)
    return F, C

def read_sequences(path):
    """sequences per locus read (column read_at). Read by hand: the read-name column of a deep sample exceeds csv's field limit."""
    sq = collections.Counter()
    with open(path) as f:
        i = f.readline().rstrip('\n').split('\t').index('read_at')
        for l in f: sq[l.split('\t', i + 1)[i]] += 1
    return sq

def selftest():
    import tempfile
    d = tempfile.mkdtemp(); p = os.path.join(d, 'x.gtf')
    open(p, 'w').write('c\tx\ttranscript\t1\t300\t.\t+\t.\tgene_id "A"; transcript_id "a1"; reads "5"; ref_gene_name "N1";\n'
                       'c\tx\texon\t1\t100\t.\t+\t.\tgene_id "A"; transcript_id "a1";\nc\tx\texon\t201\t300\t.\t+\t.\tgene_id "A"; transcript_id "a1";\n'
                       'c\tx\ttranscript\t1\t300\t.\t+\t.\tgene_id "A"; transcript_id "a2"; reads "2";\n'
                       'c\tx\texon\t11\t100\t.\t+\t.\tgene_id "A"; transcript_id "a2";\nc\tx\texon\t201\t290\t.\t+\t.\tgene_id "A"; transcript_id "a2";\n'
                       'c\tx\ttranscript\t1001\t1300\t.\t+\t.\tgene_id "B"; transcript_id "b1"; reads "3";\nc\tx\texon\t1001\t1300\t.\t+\t.\tgene_id "B"; transcript_id "b1";\n')
    g = load_gtf(p)
    assert g['A']['tx'] == 2 and len(g['A']['chains']) == 1 and g['A']['reads'] == 7 and g['A']['name'] == 'N1' and g['B']['spliced'] == 0
    cp = [dict(family_id='F0', gene_id=x, chrom='c', start='1', end='9', strand='+') for x in ('A', 'B', 'A', 'Z')]
    F, C = report(g, cp)
    assert F == [['F0', 3, 1, 3, 1, '10', 1, 'N1']], F
    assert len(C) == 3 and C[2][7:] == [0, 0, '0']
    eq = {'A': dict(placed='5', max='9', em='6.5', group='G1'), 'B': dict(placed='0', max='4', em='2.5', group='G1')}
    F, C = report(g, cp, eq)
    assert F[0][8:] == [5, 1, 1, 2] and C[0][10:] == [5, 9, '6.5', 'G1'] and C[2][10:] == [0, 0, '0.0', '.'], (F, C)
    F, C = report(g, cp, None, {'A': 3})
    assert F[0][-1] == 3 and C[0][-1] == 3 and C[1][-1] == 0, (F, C)
    sp = os.path.join(d, 's.tsv'); open(sp, 'w').write('group\tloci\tread_at\tcandidate\treads\talleles\tp\tambiguous_reads\tread_names\tisoform\n' + 'G\tA,B\tA\t1\t9\treference\tNA\t0\t' + 'r' * 300000 + '\t1\nG\tA,B\tA\t2\t3\t5SA\t0.001\t0\tx,y,z\t1\n')
    assert read_sequences(sp) == {'A': 2}
    # a locus cut into two pieces: A.p1 holds a1 / a2 (exons 1-100, 201-300), A.p2 would hold a transcript at 5001-5100 (none here)
    g = load_gtf(p); cp2 = [dict(family_id='F0', gene_id='A.p1', chrom='c', start='0', end='300', strand='+', exons='0-100,200-300'),
                            dict(family_id='F0', gene_id='A.p2', chrom='c', start='5000', end='5100', strand='+', exons='5000-5100')]
    pieces(g, cp2); F, C = report(g, cp2, None, {'A': 3})
    assert C[0][7:] == [2, 1, '7', 3] and C[1][7:] == [0, 0, '0', 0] and F[0][3] == 2, (F, C)
    print('o2_family_report selftest ok')

def main():
    a = sys.argv[1:]
    if '--selftest' in a: return selftest()
    eqp = None
    if '--eq' in a: i = a.index('--eq'); eqp = a[i + 1]; del a[i:i + 2]
    sqp = None
    if '--sequences' in a: i = a.index('--sequences'); sqp = a[i + 1]; del a[i:i + 2]
    if len(a) != 3: sys.exit(__doc__)
    gtf, cop, out = a
    sq = read_sequences(sqp) if sqp else None
    eq = {r['locus']: r for r in csv.DictReader(open(eqp), delimiter='\t')} if eqp else None
    genes = load_gtf(gtf); copies = list(csv.DictReader(open(cop), delimiter='\t')); pieces(genes, copies)
    F, C = report(genes, copies, eq, sq)
    fh = 'family_id\tcopies\tcopies_spliced\ttranscripts\tisoforms\treads\tchroms\tnames' + ('\treads_placed\ttie_groups\tmin_expressed\tmax_expressed' if eq is not None else '') + ('\ttied_sequences' if sq is not None else '')
    ch = 'family_id\tgene_id\tchrom\tstart\tend\tstrand\tname\ttranscripts\tisoforms\treads' + ('\treads_placed\treads_max\treads_em\ttie_group' if eq is not None else '') + ('\ttied_sequences' if sq is not None else '')
    for p, h, rows in ((out + '.families.tsv', fh, F), (out + '.copies.tsv', ch, C)):
        with open(p + '.tmp', 'w') as w: w.write(h + '\n'); w.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
        os.replace(p + '.tmp', p)
    multi = [r for r in F if r[1] >= 2]
    print(f'o2: {len(F)} families ({len(multi)} with >= 2 copies: {sum(r[1] for r in multi)} copies, {sum(r[2] for r in multi)} with a spliced transcript, {sum(r[4] for r in multi)} isoforms)')

if __name__ == '__main__':
    main()
