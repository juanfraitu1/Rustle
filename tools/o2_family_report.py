#!/usr/bin/env python3
"""O2 family report of one run (driver stage `o2`, the O2 stage of `all` since 2026-11-22, register 1903).

    o2_family_report.py PREFIX.gtf PREFIX.fam.copies.tsv OUT_PREFIX [--eq PREFIX.eq.loci.tsv]

What the assembly reports for each gene family of `families`: its copies (loci), the copies with a spliced transcript,
transcripts, distinct intron chains (isoforms) and reads. With --eq (stage `eqclasses`, or `multicopy`'s
PREFIX.multicopy.eq.loci.tsv) also the reads placed on one copy alone, the tie groups and the range of copies the reads
can show expressed (copies with >= 2 placed reads .. copies with any read). No read is assigned to a copy.
Writes OUT_PREFIX.families.tsv (one row per family) and OUT_PREFIX.copies.tsv (one row per copy).
The same question for several tools, with planted truth: bench/o2_families/ (docs/O2_FAMILY_BENCHMARK.md).
--selftest runs the pure parts."""
import sys, re, csv, collections

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
    out = collections.defaultdict(lambda: {'tx': 0, 'spliced': 0, 'chains': set(), 'reads': 0.0, 'name': '.'})
    for t, e in ex.items():
        g, r, n = meta[t]; e.sort(); d = out[g]; d['tx'] += 1; d['reads'] += r
        if n != '.': d['name'] = n
        if len(e) > 1: d['spliced'] += 1; d['chains'].add(tuple((a[1], b[0]) for a, b in zip(e, e[1:])))
    return out

def report(genes, copies, eq=None):
    """copies: rows of fam.copies.tsv (dicts). -> family rows, copy rows"""
    fam = collections.defaultdict(list); seen = set()
    for r in copies:
        k = (r['family_id'], r['gene_id'])
        if k not in seen: seen.add(k); fam[r['family_id']].append(r)
    F, C = [], []
    for f in sorted(fam, key=lambda f: (-len(fam[f]), f)):
        tx = sp = iso = 0; reads = 0.0; chroms = set(); names = set(); placed = mx = 0; groups = set(); lo = hi = 0
        for r in fam[f]:
            g = genes.get(r['gene_id'], {'tx': 0, 'spliced': 0, 'chains': set(), 'reads': 0.0, 'name': '.'})
            e = (eq or {}).get(r['gene_id'])
            tx += g['tx']; sp += g['spliced'] > 0; iso += len(g['chains']); reads += g['reads']; chroms.add(r['chrom'])
            if g['name'] != '.': names.add(g['name'])
            row = [f, r['gene_id'], r['chrom'], r.get('locus_start', r['start']), r.get('locus_end', r['end']), r['strand'], g['name'], g['tx'], len(g['chains']), f"{g['reads']:g}"]
            if eq is not None:
                p, m = (int(e['placed']), int(e['max'])) if e else (0, 0)
                placed += p; mx += m; lo += p >= 2; hi += m > 0
                if e and e['group'] != '.': groups.add(e['group'])
                row += [p, m, e['em'] if e else '0.0', e['group'] if e else '.']
            C.append(row)
        row = [f, len(fam[f]), sp, tx, iso, f'{reads:g}', len(chroms), ','.join(sorted(names)) or '.']
        if eq is not None: row += [placed, len(groups), max(1, lo) if hi else 0, hi]
        F.append(row)
    return F, C

def selftest():
    import tempfile, os
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
    print('o2_family_report selftest ok')

def main():
    a = sys.argv[1:]
    if '--selftest' in a: return selftest()
    eqp = None
    if '--eq' in a: i = a.index('--eq'); eqp = a[i + 1]; del a[i:i + 2]
    gtf, cop, out = a
    eq = {r['locus']: r for r in csv.DictReader(open(eqp), delimiter='\t')} if eqp else None
    F, C = report(load_gtf(gtf), list(csv.DictReader(open(cop), delimiter='\t')), eq)
    fh = 'family_id\tcopies\tcopies_spliced\ttranscripts\tisoforms\treads\tchroms\tnames' + ('\treads_placed\ttie_groups\tmin_expressed\tmax_expressed' if eq is not None else '')
    ch = 'family_id\tgene_id\tchrom\tstart\tend\tstrand\tname\ttranscripts\tisoforms\treads' + ('\treads_placed\treads_max\treads_em\ttie_group' if eq is not None else '')
    for p, h, rows in ((out + '.families.tsv', fh, F), (out + '.copies.tsv', ch, C)):
        with open(p + '.tmp', 'w') as w: w.write(h + '\n'); w.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
        __import__('os').replace(p + '.tmp', p)
    multi = [r for r in F if r[1] >= 2]
    print(f'o2: {len(F)} families ({len(multi)} with >= 2 copies: {sum(r[1] for r in multi)} copies, {sum(r[2] for r in multi)} with a spliced transcript, {sum(r[4] for r in multi)} isoforms)')

if __name__ == '__main__':
    main()
