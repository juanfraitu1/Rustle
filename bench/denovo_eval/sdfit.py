#!/usr/bin/env python3
"""Row 1755 (docs/PREREG_sd_fit_2026-11-15.md): do SDs fit multi-copy families beyond sequence similarity, and does exon
correspondence through the SEDEF CIGAR fit genes better than the span rule. Annotation level (CAT genes, Compara, CHM13 SEDEF).

  sdfit.py [SEDEF_BED] [GENOME.fai]     sdfit.py --selftest      (env as famtruth: DENOVO_CAT_GFF, DENOVO_CAT_TSV, DENOVO_PARALOGS)
"""
import sys, os, re, collections, bisect, random, statistics as st

YOUNG = {'Homo sapiens', 'Homininae', 'Hominoidea', 'Catarrhini', 'Simiiformes'}


def merged(bl):
    out = []
    for s, e in sorted(bl):
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out


def ov(bl, s, e):
    return sum(max(0, min(y, e) - max(x, s)) for x, y in bl)


def holds(bl, s, e):
    """the span rule of --sd-damp: >= half the exon bases and >= half the span inside [s, e)"""
    tot = sum(y - x for x, y in bl); lo, hi = bl[0][0], bl[-1][1]
    return 2 * ov(bl, s, e) >= tot and 2 * max(0, min(hi, e) - max(lo, s)) >= hi - lo


def aligned_blocks(s1, e1, s2, e2, minus, cigar, prefer='D'):
    """[(a_start, a_end, b_start)] per M run, b_start = forward coordinate of the b base facing a_start (walks down if minus)"""
    ops = [(int(n), o) for n, o in re.findall(r'(\d+)([MID=X])', cigar)]
    m = sum(n for n, o in ops if o in 'M=X'); i = sum(n for n, o in ops if o == 'I'); d = sum(n for n, o in ops if o == 'D')
    fits = [g for g, x, y in (('D', d, i), ('I', i, d)) if m + x == e1 - s1 and m + y == e2 - s2]
    if not fits: return None   # CIGAR does not span the stated intervals
    a_gap = prefer if prefer in fits else fits[0]   # both fit when insert and delete totals are equal: the file's convention decides
    out, a, b = [], s1, 0
    for n, o in ops:
        if o in 'M=X':
            out.append((a, a + n, e2 - 1 - b if minus else s2 + b)); a += n; b += n
        elif o == a_gap: a += n
        else: b += n
    return out


def corr_bases(blocks, minus, ex_a, ex_b):
    """exon bases of a that project onto exon bases of b"""
    tot = 0
    for a0, a1, b0 in blocks:
        for x, y in ex_a:
            lo, hi = max(x, a0), min(y, a1)
            if lo >= hi: continue
            p, q = (b0 - (hi - 1 - a0), b0 - (lo - a0) + 1) if minus else (b0 + lo - a0, b0 + hi - a0)
            tot += ov(ex_b, p, q)
    return tot


def selftest():
    assert merged([(5, 9), (1, 6), (20, 30)]) == [[1, 9], [20, 30]] and holds([[10, 20]], 14, 40) and not holds([[10, 20]], 16, 40)
    b = aligned_blocks(100, 120, 1000, 1025, False, '10M5I10M')   # 5 extra bases on side 2
    assert b == [(100, 110, 1000), (110, 120, 1015)]
    assert corr_bases(b, False, [[105, 115]], [[1000, 1012], [1015, 1030]]) == 10
    assert corr_bases(b, False, [[105, 115]], [[1010, 1015]]) == 0
    r = aligned_blocks(100, 120, 1000, 1020, True, '20M')
    assert r == [(100, 120, 1019)] and corr_bases(r, True, [[100, 105]], [[1015, 1020]]) == 5 and corr_bases(r, True, [[100, 105]], [[1000, 1015]]) == 0
    assert aligned_blocks(0, 10, 0, 10, False, '7M') is None
    print('selftest ok')


def main(sedef, fai):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import famtruth as F
    G = F.load_cat(); pc = sorted(g for g, d in G.items() if d['biotype'] == 'protein_coding')
    EX = {g: merged(G[g]['exons']) for g in pc}; NB = {g: sum(y - x for x, y in EX[g]) for g in pc}
    genes_by = collections.defaultdict(list)
    for g in pc: genes_by[G[g]['chrom']].append((EX[g][0][0], EX[g][-1][1], g))
    for v in genes_by.values(): v.sort()
    gst = {c: [x[0] for x in v] for c, v in genes_by.items()}; gmx = {c: max(x[1] - x[0] for x in v) for c, v in genes_by.items()}

    def genes_in(c, s, e):
        v = genes_by.get(c, []); out = []
        for i in range(bisect.bisect_left(gst.get(c, []), s - gmx.get(c, 0)), len(v)):
            if v[i][0] >= e: break
            if v[i][1] > s and ov(EX[v[i][2]], s, e) > 0: out.append(v[i][2])
        return out

    # truth
    by_ens = collections.defaultdict(set)
    for g in pc:
        for x in G[g]['ens']: by_ens[x].add(g)
    young, anyage = set(), set()
    for line in F._open(F.PARALOGS):
        f = line.rstrip('\n').split('\t')
        if len(f) < 4 or not f[1] or f[0].startswith('Gene'): continue
        for x in by_ens.get(f[0].split('.')[0], ()):
            for y in by_ens.get(f[1].split('.')[0], ()):
                if x != y and G[x]['name'] != G[y]['name']:
                    p = tuple(sorted((x, y))); anyage.add(p)
                    if f[3] in YOUNG: young.add(p)
    # one pass over SEDEF: span links (+ flank), exon-overlap candidates (+ CORR)
    conv = collections.Counter()
    for l in open(sedef):
        f = l.rstrip('\n').split('\t')
        if l.startswith('#') or len(f) < 33: continue
        o = {k: sum(int(n) for n, c in re.findall(r'(\d+)([MID=X])', f[32]) if c in k) for k in ('M=X', 'I', 'D')}
        if o['I'] != o['D']: conv['D' if o['M=X'] + o['D'] == int(f[2]) - int(f[1]) else 'I'] += 1
        if sum(conv.values()) >= 2000: break
    prefer = conv.most_common(1)[0][0]; bystrand = collections.Counter()
    print(f'CIGAR convention over the first 2000 unambiguous rows: {dict(conv)}')
    span_flank, corr, partial, touched, badcigar, sds = {}, {}, set(), collections.defaultdict(set), 0, []
    for n, l in enumerate(open(sedef)):
        f = l.rstrip('\n').split('\t')
        if l.startswith('#') or len(f) < 33: continue
        c1, s1, e1, c2, s2, e2 = f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]); minus = f[8] != f[9]
        sds.append((c1, s1, e1)); sds.append((c2, s2, e2))
        ga, gb = genes_in(c1, s1, e1), genes_in(c2, s2, e2)
        for g in ga: touched[g].add((n, 0))
        for g in gb: touched[g].add((n, 1))
        if not ga or not gb: continue
        blocks = None
        for a in ga:
            for b in gb:
                if a == b or G[a]['name'] == G[b]['name']: continue
                p = tuple(sorted((a, b))); partial.add(p)
                if holds(EX[a], s1, e1) and holds(EX[b], s2, e2):
                    sp = lambda g, s, e: (e - s) - max(0, min(EX[g][-1][1], e) - max(EX[g][0][0], s))
                    span_flank[p] = max(span_flank.get(p, 0), min(sp(a, s1, e1), sp(b, s2, e2)))
                if blocks is None:
                    blocks = aligned_blocks(s1, e1, s2, e2, minus, f[32], prefer) or []
                    badcigar += not blocks
                c = corr_bases(blocks, minus, EX[a], EX[b]) / min(NB[a], NB[b]); corr[p] = max(corr.get(p, 0.0), c)
                if p in young and holds(EX[a], s1, e1) and holds(EX[b], s2, e2): bystrand[(minus, c >= 0.5)] += 1
    print(f'genes {len(pc)}; young truth pairs {len(young)}; any-age pairs {len(anyage)}; SD pairs with unusable CIGAR among those walked: {badcigar}')

    # A1 flank
    fl = sorted(span_flank[p] for p in span_flank if p in young)
    q = lambda v, x: v[min(len(v) - 1, int(x * len(v)))]
    print(f'A1 flank of span-linked young pairs (n {len(fl)}): q10 {q(fl, .1)} q25 {q(fl, .25)} median {st.median(fl):.0f} q75 {q(fl, .75)}; '
          f'>= 5 kb {sum(x >= 5000 for x in fl) / len(fl):.3f}  ->  {"PASS" if st.median(fl) >= 5000 else "FAIL"} (bar median >= 5000)')

    f2 = sorted(span_flank[p] for p in span_flank if p in young and corr.get(p, 0) >= 0.5)
    print(f'A1 post hoc, only pairs that are exon counterparts (CORR >= 0.5; n {len(f2)}): q10 {q(f2, .1)} q25 {q(f2, .25)} median {st.median(f2):.0f}; >= 5 kb {sum(x >= 5000 for x in f2) / len(f2):.3f}')
    print('check: span-holding (SD pair, young pair) instances with CORR >= 0.5, same strand ' + f'{bystrand[(False, True)]} of {bystrand[(False, True)] + bystrand[(False, False)]}' + ', opposite strand ' + f'{bystrand[(True, True)]} of {bystrand[(True, True)] + bystrand[(True, False)]}')
    # A2 circular-shift null, gene level
    ygenes = {g for p in young for g in p}; clen = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(fai)}

    def share(iv):
        by = collections.defaultdict(list)
        for c, s, e in iv: by[c].append((s, e))
        ins = set()
        for c, v in by.items():
            v.sort(); stt = [x[0] for x in v]; mx = max(e - s for s, e in v)
            for lo, hi, g in genes_by.get(c, []):
                for i in range(bisect.bisect_left(stt, lo - mx), len(v)):
                    if v[i][0] >= hi: break
                    if v[i][1] > lo and holds(EX[g], *v[i]): ins.add(g); break
        return len(ins & ygenes) / len(ygenes), len(ins - ygenes) / (len(pc) - len(ygenes))
    obs = share(sds); random.seed(1); null = []
    for _ in range(100):
        off = {c: random.randrange(L) for c, L in clen.items()}
        null.append(share([(c, (s + off[c]) % (clen[c] - (e - s)), (s + off[c]) % (clen[c] - (e - s)) + e - s) for c, s, e in sds if c in clen and e - s < clen[c]]))
    ny = [x[0] for x in null]
    print(f'A2 genes with a young paralog inside an SD: observed {obs[0]:.3f} (other genes {obs[1]:.3f}); shifted SDs mean {st.mean(ny):.3f} max {max(ny):.3f} '
          f'(other genes mean {st.mean(x[1] for x in null):.3f}); empirical p < {1 / 101:.3f}' if max(ny) < obs[0] else f'A2 not separated: observed {obs[0]:.3f} null max {max(ny):.3f}')

    # B residual
    cls = collections.Counter()
    for a, b in young:
        if (a, b) in span_flank: k = 'span-linked'
        elif (a, b) in partial: k = 'partial (same SD pair, span rule fails)'
        elif touched[a] or touched[b]: k = 'one-sided / unrelated SDs'
        elif (len(EX[a]) == 1) != (len(EX[b]) == 1): k = 'no SD, retrocopy-like'
        else: k = 'no SD'
        cls[k] += 1
    print('B young truth pairs by class: ' + '; '.join(f'{k} {v} ({v / len(young):.3f})' for k, v in cls.most_common()))

    # C span rule vs exon correspondence
    def line(name, picked):
        t = len(picked & young)
        print(f'  {name}: pairs {len(picked)} young-true {t} precision {t / max(1, len(picked)):.3f} reach {t / len(young):.3f} | any-age precision {len(picked & anyage) / max(1, len(picked)):.3f}')
        return t / max(1, len(picked)), t
    print('C linked protein-coding pairs')
    ps, ts = line('span rule', set(span_flank))
    for thr in (0.2, 0.5, 0.8):
        r = line(f'CORR >= {thr}', {p for p, v in corr.items() if v >= thr})
        if thr == 0.5: pc_, tc = r
    line('span rule and CORR >= 0.5', {p for p in span_flank if corr.get(p, 0) >= 0.5})
    print(f'C verdict: precision {pc_:.3f} vs {ps:.3f} + 0.10, young-true {tc} vs {ts}  ->  {"PASS" if pc_ >= ps + 0.10 and tc >= ts else "FAIL"}')
    cv = sorted(corr[p] for p in span_flank if p in young); cf = sorted(corr[p] for p in span_flank if p not in young)
    print(f'  CORR of span-linked pairs: young-true median {st.median(cv):.2f} (q10 {q(cv, .1):.2f}); not young-true median {st.median(cf):.2f} (q90 {q(cf, .9):.2f})')

    # post hoc (not registered): each gene's counterpart, and families as connected components of the CORR links
    yg = {g for p in young for g in p}
    for thr in (0.5, 0.8):
        ed = [p for p, v in corr.items() if v >= thr]
        has = {g for p in ed if p in young for g in p}
        par = {}
        def find(x):
            while par.setdefault(x, x) != x: par[x] = par[par[x]]; x = par[x]
            return x
        for a, b in ed: par[find(a)] = find(b)
        comp = collections.defaultdict(list)
        for g in list(par): comp[find(g)].append(g)
        cp = {tuple(sorted((a, b))) for v in comp.values() for i, a in enumerate(v) for b in v[i + 1:]}
        t = len(cp & young)
        print(f'post hoc CORR >= {thr}: young-paralog genes with an exon-corresponding young partner {len(has)} of {len(yg)} ({len(has) / len(yg):.3f}); '
              f'components {len(comp)} largest {max(map(len, comp.values()))}: pairs {len(cp)} young-true {t} precision {t / len(cp):.3f} recall {t / len(young):.3f} any-age precision {len(cp & anyage) / len(cp):.3f}')


if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    else: main(sys.argv[1] if len(sys.argv) > 1 else '/anvme/workspace/v138da13-rustle/sedef_dl/incoming/sedef_send/final_human.bed',
               sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser('~/data/chm13/chm13v2.0.fa.fai'))
