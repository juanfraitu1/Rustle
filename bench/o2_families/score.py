#!/usr/bin/env python3
"""Score assemblies of planted family reads against the planted truth (O2 family benchmark).

    score.py PLANT_DIR SCENARIO LABEL=GTF [LABEL=GTF ...] [--out TABLE.tsv] [--families]

A transcript belongs to the planted copy it shares the most exon bases with (same strand, or either strand when it has
none; >= 50 bases); otherwise it is `elsewhere`. Per tool:
  copies      expressed copies with a transcript (recall = / expressed copies); copies with a transcript that are
              silent (planted 0 reads) are false copies: precision = expressed copies called / copies called
  isoforms    planted isoforms (>= 1 read) whose intron chain a transcript of that copy reproduces exactly (a one-exon
              isoform: a one-exon transcript covering >= 80% of it and inside it by >= 80%); recall and, over the tool's
              transcripts on planted copies, precision
  joined      distinct transcripts with >= 50 exon bases on two or more copies
  elsewhere   transcripts on no planted copy
--families adds one row per family. --selftest runs the pure parts."""
import sys, os, re, gzip, collections

def ov(a, b):
    """bases shared by two sorted lists of half-open intervals"""
    n = i = j = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo: n += hi - lo
        if a[i][1] < b[j][1]: i += 1
        else: j += 1
    return n

def chain(ex):
    return tuple((a[1], b[0]) for a, b in zip(ex, ex[1:]))

def load_gtf(path):
    """-> {tid: (chrom, strand, sorted half-open exons, gene_id)}"""
    tx = collections.defaultdict(list); meta = {}
    op = gzip.open if path.endswith('.gz') else open
    for l in op(path, 'rt'):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon' or int(f[4]) < int(f[3]): continue
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if not t: continue
        g = re.search(r'gene_id "([^"]+)"', f[8])
        tx[t.group(1)].append((int(f[3]) - 1, int(f[4]))); meta[t.group(1)] = (f[0], f[6], g.group(1) if g else '.')
    return {t: (meta[t][0], meta[t][1], sorted(e), meta[t][2]) for t, e in tx.items()}

def score(truth, planted, fam, tool):
    """truth {tid: (chrom, strand, exons, copy)}; planted {tid: reads}; fam {copy: family}; tool as load_gtf.
    -> {family or 'ALL': Counter}"""
    cex = collections.defaultdict(list); cloc = {}
    for t, (c, st, ex, g) in truth.items(): cex[g] += ex; cloc[g] = (c, st)
    for g in cex: cex[g] = merge(cex[g])
    by = collections.defaultdict(list)
    for g, (c, st) in cloc.items(): by[c].append(g)
    creads = collections.Counter()
    for t, n in planted.items(): creads[truth[t][3]] += n
    iso = collections.defaultdict(list)
    for t, n in planted.items():
        if n > 0: iso[truth[t][3]].append(truth[t][2])
    R = collections.defaultdict(collections.Counter); called = set(); got = set()
    def add(g, k, n=1):
        R['ALL'][k] += n
        if g is not None: R[fam[g]][k] += n
    for t, (c, st, ex, _) in tool.items():
        hits = {}
        for g in by.get(c, ()):
            if st in '+-' and st != cloc[g][1]: continue
            if cex[g][-1][1] <= ex[0][0] or cex[g][0][0] >= ex[-1][1]: continue
            n = ov(ex, cex[g])
            if n >= 50: hits[g] = n
        if not hits: add(None, 'elsewhere'); continue
        g = max(hits, key=lambda k: (hits[k], k)); called.add(g); add(g, 'tx_on_copies')
        if len(hits) > 1: add(g, 'joined')
        ok = False
        for i, pex in enumerate(iso.get(g, ())):
            if len(pex) > 1: m = len(ex) > 1 and chain(ex) == chain(pex)
            else:
                o = ov(ex, pex); m = len(ex) == 1 and o >= 0.8 * (pex[0][1] - pex[0][0]) and o >= 0.8 * (ex[0][1] - ex[0][0])
            if m: ok = True; got.add((g, i))
        if ok: add(g, 'tx_exact')
    for g in cex:
        e = creads[g] > 0
        add(g, 'expressed' if e else 'silent')
        if g in called: add(g, 'expressed_called' if e else 'silent_called')
        add(g, 'isoforms', len(iso.get(g, ())))
        add(g, 'isoforms_found', sum(1 for i in range(len(iso.get(g, ()))) if (g, i) in got))
    R['_called'] = called
    return R

def merge(iv):
    o = []
    for s, e in sorted(iv):
        if o and s <= o[-1][1]: o[-1] = (o[-1][0], max(o[-1][1], e))
        else: o.append((s, e))
    return o

COLS = 'scenario tool family expressed copies_found copy_recall silent false_copies copy_precision isoforms isoforms_found isoform_recall tx_on_copies tx_exact isoform_precision joined elsewhere'.split()
def row(sc, lab, f, c):
    r = lambda a, b: f'{a / b:.3f}' if b else 'NA'
    called = c['expressed_called'] + c['silent_called']
    return [sc, lab, f, c['expressed'], c['expressed_called'], r(c['expressed_called'], c['expressed']), c['silent'], c['silent_called'], r(c['expressed_called'], called),
            c['isoforms'], c['isoforms_found'], r(c['isoforms_found'], c['isoforms']), c['tx_on_copies'], c['tx_exact'], r(c['tx_exact'], c['tx_on_copies']), c['joined'], c['elsewhere']]

def selftest():
    assert ov([(0, 10), (20, 30)], [(5, 25)]) == 10 and merge([(5, 9), (0, 6), (20, 30)]) == [(0, 9), (20, 30)]
    truth = {'a1': ('c', '+', [(0, 100), (200, 300)], 'A'), 'b1': ('c', '+', [(1000, 1100), (1200, 1300)], 'B'), 'm1': ('c', '-', [(5000, 5400)], 'M')}
    planted = {'a1': 10, 'b1': 0, 'm1': 5}; fam = {'A': 'F', 'B': 'F', 'M': 'G'}
    tool = {'x': ('c', '+', [(10, 100), (200, 290)], '.'),      # exact at A
            'y': ('c', '+', [(1000, 1100), (1200, 1300)], '.'),  # at silent B: false copy, not exact
            'z': ('c', '+', [(0, 100), (1200, 1300)], '.'),      # joins A and B
            'w': ('c', '.', [(5010, 5390)], '.'),                # one exon at M, strandless
            'q': ('d', '+', [(0, 500)], '.')}                    # elsewhere
    R = score(truth, planted, fam, tool); a = R['ALL']
    assert (a['expressed'], a['expressed_called'], a['silent'], a['silent_called']) == (2, 2, 1, 1), a
    assert (a['isoforms'], a['isoforms_found'], a['tx_on_copies'], a['tx_exact'], a['joined'], a['elsewhere']) == (2, 2, 4, 2, 1, 1), a
    assert R['G']['tx_exact'] == 1 and R['F']['silent_called'] == 1 and R['_called'] == {'A', 'B', 'M'}
    assert row('S', 't', 'ALL', a)[5] == '1.000' and row('S', 't', 'ALL', a)[8] == '0.667'
    print('score selftest ok')

def load_plan(d, sc):
    """-> truth {tid: (chrom, strand, exons, copy)} of the copies in scenario sc, planted {tid: reads}, fam {copy: family}"""
    T = load_gtf(f'{d}/loci.gtf'); planted, fam = {}, {}; hdr = None
    for l in open(f'{d}/plan.tsv'):
        f = l.rstrip('\n').split('\t')
        if hdr is None: hdr = f; continue
        fam[f[1]] = f[0]
        if f[hdr.index(sc)] != 'NA': planted[f[3]] = int(f[hdr.index(sc)])
    keep = {T[t][3] for t in planted}; truth = {t: v for t, v in T.items() if v[3] in keep}
    return truth, {t: planted.get(t, 0) for t in truth}, fam

def main():
    a = sys.argv[1:]
    if '--selftest' in a: return selftest()
    out = None; fams = '--families' in a
    if fams: a.remove('--families')
    if '--out' in a: i = a.index('--out'); out = a[i + 1]; del a[i:i + 2]
    d, sc = a[0], a[1]
    T = load_gtf(f'{d}/loci.gtf'); truth = {t: (c, st, ex, g) for t, (c, st, ex, g) in T.items()}
    planted, fam = {}, {}; hdr = None
    for l in open(f'{d}/plan.tsv'):
        f = l.rstrip('\n').split('\t')
        if hdr is None: hdr = f; continue
        fam[f[1]] = f[0]
        if f[hdr.index(sc)] != 'NA': planted[f[3]] = int(f[hdr.index(sc)])
    keep = {truth[t][3] for t in planted}; truth = {t: v for t, v in truth.items() if v[3] in keep}
    planted = {t: planted.get(t, 0) for t in truth}
    rows = []
    for x in a[2:]:
        lab, path = x.split('=', 1)
        if not os.path.exists(path): rows.append([sc, lab, 'MISSING'] + [''] * (len(COLS) - 3)); continue
        R = score(truth, planted, fam, load_gtf(path))
        rows.append(row(sc, lab, 'ALL', R['ALL']))
        if fams: rows += [row(sc, lab, f, R[f]) for f in sorted(R) if f not in ('ALL', '_called')]
    txt = '\t'.join(COLS) + '\n' + ''.join('\t'.join(map(str, r)) + '\n' for r in rows)
    sys.stdout.write(txt)
    if out: open(out, 'w').write(txt)

if __name__ == '__main__':
    main()
