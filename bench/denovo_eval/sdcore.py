#!/usr/bin/env python3
"""Row 1756 (docs/PREREG_sd_core_label_2026-11-15.md): opt-in label on pairs of copies inside one family.
`sd_core` = one SEDEF pair holds both copies by the --sd-damp span rule AND projects, through its CIGAR, exon bases of one onto
exon bases of the other for >= 0.5 of the smaller copy's exon bases. A table only: it changes no family.

  sdcore.py COPIES_TSV SEDEF_BED OUT.tsv     one row per family copy pair that has exon bases in the two copies of one SEDEF pair
  sdcore.py score COPIES_TSV SDCORE.tsv [famtruth species flags] --gtf RUN.families.gtf     precision of flagged vs other family pairs
  sdcore.py --selftest
"""
import sys, os, re, collections, bisect, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sdfit import merged, ov, holds, aligned_blocks

CORR_MIN = 0.5
CIGAR = re.compile(r'(\d+[MID=X])+$')


def project(blocks, starts, minus, ex):
    """the intervals on the other copy that the exon blocks ex align to (as sdfit.corr_bases, one bisect per exon)"""
    out = []
    for x, y in ex:
        i = max(0, bisect.bisect_right(starts, x) - 1)
        while i < len(blocks) and blocks[i][0] < y:
            a0, a1, b0 = blocks[i]; lo, hi = max(x, a0), min(y, a1)
            if lo < hi: out.append((b0 - (hi - 1 - a0), b0 - (lo - a0) + 1) if minus else (b0 + lo - a0, b0 + hi - a0))
            i += 1
    return out


def load_copies(path):
    """(family, copy_idx) -> (tid, chrom, merged exon blocks)"""
    out = {}
    with open(path) as f:
        h = {k: i for i, k in enumerate(f.readline().rstrip('\n').split('\t'))}
        for l in f:
            p = l.rstrip('\n').split('\t')
            ex = merged(tuple(map(int, b.split('-'))) for b in p[h['exons']].split(',') if b)
            if ex: out[(p[h['family_id']], p[h['copy_idx']])] = (p[h['tid']], p[h['chrom']], ex)
    return out


def label(copies, sedef_lines):
    """{(key_a, key_b): [span_held, corr, sd_identity]} best over SEDEF pairs, for same-family copy pairs"""
    by = collections.defaultdict(list)
    for k, (_, c, ex) in copies.items(): by[c].append((ex[0][0], ex[-1][1], k))
    for v in by.values(): v.sort()
    st = {c: [x[0] for x in v] for c, v in by.items()}; mx = {c: max(x[1] - x[0] for x in v) for c, v in by.items()}

    def inside(c, s, e):
        v = by.get(c, []); out = []
        for i in range(bisect.bisect_left(st.get(c, []), s - mx.get(c, 0)), len(v)):
            if v[i][0] >= e: break
            if v[i][1] > s and ov(copies[v[i][2]][2], s, e) > 0: out.append(v[i][2])
        return out
    res, cig, idc, prefer, bad = {}, None, None, None, 0
    for l in sedef_lines:
        f = l.rstrip('\n').split('\t')
        if l.startswith('#') or len(f) < 33 or not f[1].isdigit(): continue
        if cig is None:  # per file, from its first row: CIGAR is col 33 (CHM13, gorilla) or 34 (HG002 build); fracMatch col 21 or 22
            cig = next((i for i in (32, 33) if i < len(f) and CIGAR.match(f[i])), None); idc = 20 if '.' in f[20] else 21
            if cig is None: sys.exit('sdcore: this SEDEF BED has no CIGAR column (33 or 34); the UCSC track does not carry one, use a native SEDEF final.bed')
        c1, s1, e1, c2, s2, e2 = f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]); minus = f[8] != f[9]
        ka, kb = inside(c1, s1, e1), inside(c2, s2, e2)
        fb = {}
        for k in kb: fb.setdefault(k[0], []).append(k)
        blocks = None
        for a in ka:
            proj = None
            for b in fb.get(a[0], ()):
                if a == b: continue
                if blocks is None:
                    if prefer is None:  # the file's insert/delete convention, from the first row where the two totals differ
                        o = {k: sum(int(n) for n, c in re.findall(r'(\d+)([MID=X])', f[cig]) if c in k) for k in ('M=X', 'I', 'D')}
                        if o['I'] != o['D']: prefer = 'D' if o['M=X'] + o['D'] == e1 - s1 else 'I'
                    blocks = aligned_blocks(s1, e1, s2, e2, minus, f[cig], prefer or 'D') or []
                    bad += not blocks; starts = [x[0] for x in blocks]
                ea, eb = copies[a][2], copies[b][2]
                if proj is None: proj = project(blocks, starts, minus, ea)
                corr = sum(ov(eb, p, q) for p, q in proj) / min(sum(y - x for x, y in ea), sum(y - x for x, y in eb))
                held = holds(ea, s1, e1) and holds(eb, s2, e2)
                key = (a, b) if a < b else (b, a); r = res.setdefault(key, [0, 0.0, 0.0, 0])
                r[0] |= held; r[1] = max(r[1], corr)
                if held and corr >= CORR_MIN: r[3] = 1; r[2] = max(r[2], float(f[idc]))
    return res, bad


def run(copies_path, sedef, out):
    copies = load_copies(copies_path); res, bad = label(copies, open(sedef))
    with open(out, 'w') as o:
        o.write('family_id\tcopy_a\tcopy_b\ttid_a\ttid_b\tsd_span_held\tsd_exon_corr\tsd_core\tsd_core_identity\n')
        for (a, b), (held, corr, idn, core) in sorted(res.items()):
            o.write(f'{a[0]}\t{a[1]}\t{b[1]}\t{copies[a][0]}\t{copies[b][0]}\t{int(held)}\t{corr:.3f}\t{core}\t{idn if core else "NA"}\n')
    n = sum(r[3] for r in res.values())
    print(f'sd_core pairs {n} in {len({k[0][0] for k, r in res.items() if r[3]})} families, of {len(res)} family copy pairs sharing a SEDEF pair (unusable CIGARs {bad})')


def score(a):
    """flagged vs other within-family protein-coding gene pairs, judged as pairscore judges"""
    import csv, famtruth as F
    copies_path, table = a[0], a[1]; gtf = a[a.index('--gtf') + 1]
    G = F.genes_from_argv(a, a[a.index('--species') + 1] if '--species' in a else 'human'); idx = F.exon_index(G)
    expressed = F.expressed_from_gtf(gtf, idx); truth, universe = F.truth_pairs(G, F.need(F.PARALOGS, 'DENOVO_PARALOGS'), F.levels('simiiformes'))
    truth = {tuple(sorted(p)) for p in truth if p[0] in expressed and p[1] in expressed}
    gene, fams = {}, collections.defaultdict(set)
    for r in csv.DictReader(open(copies_path), delimiter='\t'):
        g = F.map_copy(idx, r['chrom'], r['strand'] if r['strand'] in '+-' else '+', F.blocks(r['exons']))
        if g: gene[(r['family_id'], r['copy_idx'])] = g; fams[r['family_id']].add(g)
    ok = lambda x, y: (G[x]['biotype'] == 'protein_coding' and G[y]['biotype'] == 'protein_coding' and x in universe and y in universe
                       and G[x]['name'] != G[y]['name'])
    pred = {(x, y) for v in fams.values() for x, y in itertools.combinations(sorted(v), 2) if ok(x, y)}
    flag = set()
    for r in csv.DictReader(open(table), delimiter='\t'):
        x, y = gene.get((r['family_id'], r['copy_a'])), gene.get((r['family_id'], r['copy_b']))
        if r['sd_core'] == '1' and x and y and x != y: flag.add(tuple(sorted((x, y))))
    flag &= pred; rest = pred - flag
    pr = lambda s: len(s & truth) / len(s) if s else float('nan')
    anyp = set()   # Compara pairs of any age, as pairscore's `related` (without the HGNC groups)
    for r in csv.reader(F._open(F.PARALOGS), delimiter='\t'):
        if len(r) >= 2 and r[1] and not r[0].startswith('Gene'): anyp.add(frozenset((r[0].split('.')[0], r[1].split('.')[0])))
    rel = lambda s: sum(1 for x, y in s if (x, y) in truth or any(frozenset((p, q)) in anyp for p in G[x]['ens'] for q in G[y]['ens'])) / len(s) if s else float('nan')
    tf, tr = len(flag & truth), len(rest & truth)
    print('run\tfamily_pairs\tprecision_all\tflagged\tflagged_true\tprecision_flagged\tother\tother_true\tprecision_other\tflagged_share_of_true\tflagged_true_or_older_relative\tother_true_or_older_relative')
    print(f'{copies_path}\t{len(pred)}\t{pr(pred):.3f}\t{len(flag)}\t{tf}\t{pr(flag):.3f}\t{len(rest)}\t{tr}\t{pr(rest):.3f}\t{tf / max(1, tf + tr):.3f}\t{rel(flag):.3f}\t{rel(rest):.3f}')


def selftest():
    cop = {('F', '0'): ('a', 'c1', [[100, 200], [300, 400]]), ('F', '1'): ('b', 'c2', [[1100, 1200], [1300, 1400]]),
           ('F', '2'): ('c', 'c2', [[1500, 1600]]), ('G', '0'): ('d', 'c2', [[1100, 1200]])}
    row = lambda s1, e1, s2, e2, cg: '\t'.join(['c1', str(s1), str(e1), 'c2', str(s2), str(e2), 'S', '1', '+', '+'] + ['x'] * 10 + ['0.99'] + ['x'] * 11 + [cg, '0.9'])
    res, bad = label(cop, [row(0, 700, 1000, 1700, '700M')])
    assert bad == 0 and res[(('F', '0'), ('F', '1'))] == [1, 1.0, 0.99, 1], res      # counterpart: held and exon-matched
    assert res[(('F', '0'), ('F', '2'))][1:] == [0.0, 0.0, 0] and (('F', '0'), ('G', '0')) not in res  # same SD, other position; other family
    res, _ = label(cop, [row(0, 700, 1040, 1740, '700M')])                             # shifted 40 bp: 60 of 100 exon bases match
    assert res[(('F', '0'), ('F', '1'))][1] == 0.6
    print('selftest ok')


if __name__ == '__main__':
    a = sys.argv[1:]
    if a == ['--selftest']: selftest()
    elif a and a[0] == 'score': score(a[1:])
    elif len(a) == 3: run(*a)
    else: sys.exit(__doc__)
