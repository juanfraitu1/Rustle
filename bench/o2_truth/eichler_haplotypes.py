#!/usr/bin/env python3
"""The Eichler lab's 10-point filter (Guitart 2024) on an individual's own haplotypes (register 1972).

    eichler_haplotypes.py genes CHM13.fa OUT.fa OUT.tsv            family genes' genomic sequences (env as famtruth)
    eichler_haplotypes.py score GENES.tsv MAT.paf PAT.paf MAT.bam PAT.bam CHM13.bam [MIN_LEN=1000]

score: per read its primary on each haplotype (the paper's mapping), the paralog under it (the reference copy whose
sequence has the most matching bases on the haplotype inside the read's span, from the asm5 PAFs), and the rule: same paralog on both, or one
mapping = retained; different paralogs = retained for the higher AS when it leads by >= 10, else ambiguous. Beside it
the read on CHM13: the paralog under the primary (>= 50 exon bases), near / exact tie on ms, AS lead of the primary."""
import sys, os, collections, bisect

def genes(fa, out_fa, out_tsv):
    import pysam
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval')); import famtruth as F
    G = F.load_genes('human'); pairs, _ = F.truth_pairs(G, F.PARALOGS, set(F.LEVELS['simiiformes'])); fa = pysam.FastaFile(fa); n = 0
    with open(out_fa, 'w') as w, open(out_tsv, 'w') as t:
        t.write('gene\tname\tfamily\tchrom\tstart\tend\tstrand\texons\n')
        for i, comp in enumerate(F.components(pairs)):
            for g in sorted(comp):
                v = G[g]
                if not v['exons'] or v['chrom'] not in fa.references: continue
                s, e = v['exons'][0][0], v['exons'][-1][1]; w.write(f'>{g}\n{fa.fetch(v["chrom"], s, e)}\n'); n += 1
                t.write(f"{g}\t{v['name']}\tF{i}\t{v['chrom']}\t{s}\t{e}\t{v['strand']}\t{','.join(f'{a}-{b}' for a, b in v['exons'])}\n")
    print(f'{n} genes written')

def load_paf(p):
    """-> {contig: sorted [(start, end, gene, run starts, run ends, cumulative matches)]} of the hits covering >= half of the gene
    sequence; the runs are the hit's stretches of matching bases on the haplotype (the `=` of its cigar)"""
    import re
    by = collections.defaultdict(list)
    for l in open(p):
        f = l.rstrip('\n').split('\t')
        if (int(f[3]) - int(f[2])) < 0.5 * int(f[1]): continue
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        if cg is None: continue
        t = int(f[7]); rs = []; re_ = []; cum = [0]
        for n, o in re.findall(r'(\d+)([=XIDMN])', cg):
            n = int(n)
            if o == '=': rs.append(t); re_.append(t + n); cum.append(cum[-1] + n)
            if o in '=XDMN': t += n
        by[f[5]].append((int(f[7]), int(f[8]), f[0], rs, re_, cum))
    for v in by.values(): v.sort(key=lambda x: (x[0], x[1]))
    return by, {c: max(x[1] - x[0] for x in v) for c, v in by.items()}

def window_matches(rs, re_, cum, s, e):
    """matching bases of a hit inside [s, e)"""
    i = bisect.bisect_right(re_, s); j = bisect.bisect_left(rs, e)
    if i >= j: return 0
    return cum[j] - cum[i] - max(0, s - rs[i]) - max(0, re_[j - 1] - e)

def label(hits, mx, c, s, e):
    """the reference copy with the most matching bases on the haplotype INSIDE the read's span (a hit must cover >= 80% of it)"""
    v = hits.get(c)
    if not v: return None
    k = bisect.bisect_left(v, (s - mx[c],)); best = (0, None)
    while k < len(v) and v[k][0] < e:
        a, b, g, rs, re_, cum = v[k]; k += 1
        if min(b, e) - max(a, s) >= 0.8 * (e - s):
            m = window_matches(rs, re_, cum, s, e)
            if m > best[0]: best = (m, g)
    return best[1]

def primaries(bam, minlen):
    import pysam
    d = {}
    for r in pysam.AlignmentFile(bam, check_sq=False):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        if (r.infer_read_length() or 0) < minlen: continue
        d[r.query_name] = (r.reference_name, r.reference_start, r.reference_end, r.get_tag('AS'))
    return d

def rule(m, p):
    """m, p = (paralog, AS) or None per haplotype -> ('retained', paralog) | ('ambiguous', None) | ('unmapped', None)"""
    if m is None and p is None: return 'unmapped', None
    if m is None or p is None: return 'retained', (m or p)[0]
    if m[0] == p[0]: return 'retained', m[0]
    hi, lo = (m, p) if m[1] >= p[1] else (p, m)
    return ('retained', hi[0]) if hi[1] - lo[1] >= 10 else ('ambiguous', None)

def score(gtsv, mpaf, ppaf, mbam, pbam, cbam, minlen=1000):
    import pysam
    minlen = int(minlen); G = {}; byc = collections.defaultdict(list)
    for l in list(open(gtsv))[1:]:
        g, name, fam, c, s, e, st, ex = l.rstrip('\n').split('\t'); exs = [tuple(map(int, x.split('-'))) for x in ex.split(',')]; G[g] = (name, fam); byc[c].append((int(s), int(e), exs, g))
    H = {}; 
    for h, paf, bam in (('mat', mpaf, mbam), ('pat', ppaf, pbam)):
        hits, mx = load_paf(paf); pr = primaries(bam, minlen); H[h] = {q: (label(hits, mx, c, s, e), a) for q, (c, s, e, a) in pr.items()}
        print(f'{h}: {len(pr)} reads with a primary, {sum(1 for v in H[h].values() if v[0])} on a labelled family locus')
    # the single reference: primary's paralog, ms tie status, AS lead
    best = collections.defaultdict(lambda: [-10 ** 9, -10 ** 9, None, None, [], []])    # best ms, primary AS, primary paralog, primary ms, all ms, other AS
    for r in pysam.AlignmentFile(cbam, check_sq=False):
        if r.is_unmapped or r.is_supplementary: continue
        b = best[r.query_name]; ms = r.get_tag('ms'); a = r.get_tag('AS'); b[4].append(ms)
        if r.is_secondary: b[5].append(a); continue
        b[1] = a; b[3] = ms; bl = r.get_blocks(); lab = (49, None)
        for s, e, exs, g in byc.get(r.reference_name, ()):
            if s < bl[-1][1] and e > bl[0][0]:
                o = sum(max(0, min(y, q) - max(x, p)) for x, y in bl for p, q in exs)
                if o > lab[0]: lab = (o, g)
        b[2] = lab[1]
    OUT = open(os.environ['EICHLER_READS_OUT'], 'w') if os.environ.get('EICHLER_READS_OUT') else None   # read, verdict, paralog: the truth of register 1974
    def report(title, keep):
        C = collections.Counter(); A = collections.Counter(); T = collections.Counter(); S = collections.Counter(); n = 0; X = collections.Counter()
        for q in set(H['mat']) | set(H['pat']):
            m, p = H['mat'].get(q), H['pat'].get(q)
            labs = {x[0] for x in (m, p) if x and x[0]}
            if not labs or not keep(labs): continue
            if (m and not m[0]) or (p and not p[0]): C['one haplotype primary off the families'] += 1
            m2 = m if m and m[0] else None; p2 = p if p and p[0] else None; v, par = rule(m2, p2); n += 1; C[v] += 1
            if OUT is not None and title == 'all truth families': OUT.write(f'{q}\t{v}\t{par or "."}\n')
            b = best.get(q)
            if b is None or b[3] is None: A[(v, 'not mapped on CHM13')] += 1; continue
            top = max(b[4]); near = sum(1 for x in b[4] if x >= top - 10) > 1 or b[3] < top; exact = sum(1 for x in b[4] if x == top) > 1
            lead = b[1] - max(b[5]) if b[5] else 10 ** 6
            if v == 'retained': A[(v, 'CHM13 primary on the same paralog' if b[2] == par else 'on another family paralog' if b[2] else 'on no family gene')] += 1
            tie = 'exact tie' if exact else 'near tie' if near else 'clear'
            if v == 'retained': X[(tie, 'same paralog' if b[2] == par else 'same family, another paralog' if b[2] and G[b[2]][1] == G[par][1] else 'another family or none')] += 1
            T[(v, tie)] += 1; S[(v, 'fires' if lead >= 10 else 'abstains')] += 1
        print(f'## {title}: {n} reads with a haplotype primary on a labelled family locus')
        for v in ('retained', 'ambiguous'):
            k = C[v]; print(f'  {v}: {k} ({k / max(1, n):.3f})')
            if v == 'retained': print('     ' + '; '.join(f'{x}: {A[(v, x)]} ({A[(v, x)] / max(1, k):.3f})' for x in ('CHM13 primary on the same paralog', 'on another family paralog', 'on no family gene', 'not mapped on CHM13')))
            print('     on the single reference, ms: ' + '; '.join(f'{x}: {T[(v, x)]} ({T[(v, x)] / max(1, k):.3f})' for x in ('clear', 'near tie', 'exact tie')) + ' | single-reference adaptation: ' + '; '.join(f'{x}: {S[(v, x)]} ({S[(v, x)] / max(1, k):.3f})' for x in ('fires', 'abstains')))
        for t in ('clear', 'near tie', 'exact tie'):
            k = sum(v for (tt, _), v in X.items() if tt == t); print(f'     retained and {t} on the single reference ({k}): ' + '; '.join(f'{x} {X[(t, x)]} ({X[(t, x)] / max(1, k):.3f})' for x in ('same paralog', 'same family, another paralog', 'another family or none')))
        print(f"  (reads with one haplotype primary off the families: {C['one haplotype primary off the families']})")
    report('all truth families', lambda labs: True)
    if OUT is not None: OUT.close()
    for fam in ('TBC1D3', 'NPIP'): report(fam, lambda labs, f=fam: any(G[g][0].startswith(f) for g in labs))

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        assert rule(('A', 100), ('A', 50)) == ('retained', 'A') and rule(('A', 100), ('B', 90)) == ('retained', 'A') and rule(('A', 100), ('B', 91)) == ('ambiguous', None) and rule(None, ('B', 5)) == ('retained', 'B') and rule(None, None)[0] == 'unmapped'
        h = {'c': [(0, 1000, 'A', [0, 600], [500, 1000], [0, 500, 900]), (0, 1000, 'B', [0], [1000], [0, 1000]), (5000, 6000, 'C', [5000], [6000], [0, 1000])]}; mx = {'c': 1000}
        assert window_matches([0, 600], [500, 1000], [0, 500, 900], 100, 900) == 700 and window_matches([0], [1000], [0, 1000], 100, 900) == 800 and window_matches([0, 600], [500, 1000], [0, 500, 900], 520, 580) == 0
        assert label(h, mx, 'c', 100, 900) == 'B' and label(h, mx, 'c', 900, 2000) is None and label(h, mx, 'd', 0, 10) is None; print('eichler_haplotypes selftest ok')
    elif sys.argv[1] == 'genes': genes(*sys.argv[2:])
    else: score(*sys.argv[2:])
