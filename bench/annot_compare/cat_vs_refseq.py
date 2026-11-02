#!/usr/bin/env python3
"""CHM13 CAT/Liftoff vs NCBI RefSeq (native CHM13 annotation), descriptive (2026-11-02, no decision rule):
  missing  : genes of one annotation with no same-strand exonic overlap in the other, by cause
             (catalogue = no Ensembl counterpart; projected-elsewhere = its ENSG is a CAT source gene at another locus;
              not-projected = ENSG known to HGNC but no CAT gene) x location (T2T-unique / SD90 / other)
  longer   : same-gene pairs where the CAT span is >= 1.5x and >= 10 kb longer: merged (collapsed_gene_ids),
             stretched (max intron >= 2x and >= 10 kb longer), isoform-extent (the rest); and the reverse
  merges   : CAT collapsed genes: RefSeq genes and our assembled copies (fam.copies) at the same exons
usage: cat_vs_refseq.py CAT_GENES.gff3 CAT_SLIM.gff3.gz REFSEQ.gtf.gz SD.bed.gz UNIQUE.bed HGNC.txt PARALOGS.tsv COPIES.tsv OUTPREFIX"""
import bisect, collections, gzip, re, sys

cat_genes, cat_slim, rs_gtf, sd_bed, uniq_bed, hgnc, paralogs, copies, out = sys.argv[1:10]
op = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
A = lambda k, a: (re.search(r'(?:^|;)\s*' + k + r'[= ]"?([^;"]+)', a) or [None, ''])[1]


def merge(iv):
    o = []
    for a, b in sorted(iv):
        if o and a <= o[-1][1]:
            o[-1][1] = max(o[-1][1], b)
        else:
            o.append([a, b])
    return o


class Idx:
    def __init__(self, items):  # items: dict id -> (chrom, strand, exons)
        self.by = collections.defaultdict(list)
        for k, (c, s, ex) in items.items():
            self.by[(c, s)].append((ex[0][0], ex[-1][1], k, ex))
        self.lo = {}
        self.ml = {}
        for k, v in self.by.items():
            v.sort()
            self.lo[k] = [x[0] for x in v]
            self.ml[k] = max(x[1] - x[0] for x in v)

    def hits(self, c, s, ex):
        v = self.by.get((c, s))
        if not v:
            return []
        lo, hi = ex[0][0], ex[-1][1]
        i = bisect.bisect_left(self.lo[(c, s)], lo - self.ml[(c, s)])
        o = []
        while i < len(v) and v[i][0] < hi:
            a, b, k, e2 = v[i]
            i += 1
            if b > lo and share(ex, e2):
                o.append(k)
        return o


def share(x, y):
    i = j = 0
    while i < len(x) and j < len(y):
        if x[i][0] < y[j][1] and y[j][0] < x[i][1]:
            return True
        if x[i][1] < y[j][1]:
            i += 1
        else:
            j += 1
    return False


def bed(p):
    d = collections.defaultdict(list)
    for l in op(p):
        f = l.split('\t')
        if len(f) >= 3 and not l.startswith(('#', 'track')):
            d[f[0]].append((int(f[1]), int(f[2])))
    return {c: merge(v) for c, v in d.items()}


def frac_in(b, c, lo, hi):
    s = 0
    for a, e in b.get(c, ()):
        if e > lo and a < hi:
            s += min(e, hi) - max(a, lo)
    return s / max(hi - lo, 1)


SD, UQ = bed(sd_bed), bed(uniq_bed)
loc = lambda c, lo, hi: 'T2T-unique' if frac_in(UQ, c, lo, hi) >= 0.5 else ('SD90' if frac_in(SD, c, lo, hi) >= 0.5 else 'other')

# ---- CAT gene attributes + slim exons
cattr = {}
for l in op(cat_genes):
    f = l.rstrip('\n').split('\t')
    if len(f) > 8 and f[2] == 'gene':
        a = f[8]
        cattr[A('ID', a)] = dict(name=A('gene_name', a), src=A('source_gene', a).split('.')[0], coll=A('collapsed_gene_ids', a),
                                 colln=A('collapsed_gene_names', a), alt=A('gene_alternate_contigs', a), bt=A('gene_biotype', a))
cg, ctx, cpar = {}, collections.defaultdict(list), {}
for l in op(cat_slim):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9:
        continue
    a = f[8]
    if f[2] == 'gene':
        cg[A('ID', a)] = [f[0], f[6], A('gene_name', a), A('gene_biotype', a)]
    elif f[2] == 'transcript':
        cpar[A('ID', a)] = A('Parent', a)
    elif f[2] == 'exon':
        ctx[A('Parent', a)].append((int(f[3]) - 1, int(f[4])))
cex, ctr = collections.defaultdict(list), collections.defaultdict(list)
for t, ex in ctx.items():
    g = cpar.get(t)
    if g in cg:
        cex[g] += ex
        ctr[g].append(sorted(ex))
CAT = {g: (cg[g][0], cg[g][1], merge(e)) for g, e in cex.items()}
cat_src = {v['src'] for v in cattr.values() if v['src']}

# ---- RefSeq
rg, rex, rtr = {}, collections.defaultdict(list), collections.defaultdict(lambda: collections.defaultdict(list))
for l in op(rs_gtf):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon':
        continue
    g, t = A('gene_id', f[8]), A('transcript_id', f[8])
    rg.setdefault(g, [f[0], f[6], A('gene_name', f[8]), set()])[3].add(f[1])
    rex[g].append((int(f[3]) - 1, int(f[4])))
    rtr[g][t].append((int(f[3]) - 1, int(f[4])))
RS = {g: (rg[g][0], rg[g][1], merge(e)) for g, e in rex.items()}
coding = {g for g, d in rtr.items() if any(t.startswith(('rna-NM_', 'rna-XM_')) for t in d)}

# ---- HGNC symbol -> ENSG
sym2e, e2sym = {}, {}
hd = None
for l in open(hgnc):
    f = l.rstrip('\n').split('\t')
    if hd is None:
        hd = {k: i for i, k in enumerate(f)}
        continue
    e = f[hd['ensembl_gene_id']] if len(f) > hd['ensembl_gene_id'] else ''
    if e:
        sym2e[f[hd['symbol']]] = e
        e2sym[e] = f[hd['symbol']]
par = collections.defaultdict(set)
for l in open(paralogs):
    f = l.rstrip('\n').split('\t')
    if len(f) >= 2 and f[1]:
        par[f[0]].add(f[1])

ci, ri = Idx(CAT), Idx(RS)
L = []
P = lambda s='': (L.append(s), print(s))

# ===== missing
P('## RefSeq genes with no CAT gene at their exons (same strand)')
tab = collections.Counter()
ex_rows = collections.defaultdict(list)
for g, (c, s, ex) in RS.items():
    if ci.hits(c, s, ex):
        continue
    nm = rg[g][2]
    e = sym2e.get(nm)
    cause = 'catalogue (no Ensembl counterpart)' if not e else ('projected elsewhere' if e in cat_src else 'not projected')
    where = loc(c, ex[0][0], ex[-1][1])
    tab[(cause, where, g in coding)] += 1
    if len(ex_rows[(cause, where)]) < 8 and g in coding:
        ex_rows[(cause, where)].append(nm)
P('| cause | location | coding | non-coding |')
P('|---|---|---|---|')
for cause in ('catalogue (no Ensembl counterpart)', 'projected elsewhere', 'not projected'):
    for where in ('T2T-unique', 'SD90', 'other'):
        P(f'| {cause} | {where} | {tab[(cause, where, True)]} | {tab[(cause, where, False)]} |')
P(f'RefSeq genes {len(RS)} (coding {len(coding)}); missing in CAT {sum(tab.values())} (coding {sum(v for k, v in tab.items() if k[2])})')
for k, v in ex_rows.items():
    P(f'  e.g. {k[0]} / {k[1]}: ' + ' '.join(v))

P('\n## CAT genes with no RefSeq gene at their exons')
tab2 = collections.Counter()
for g, (c, s, ex) in CAT.items():
    if ri.hits(c, s, ex):
        continue
    bt = cg[g][3] or 'NA'
    bt = bt if bt in ('protein_coding', 'lncRNA') else 'other'
    tab2[(bt, loc(c, ex[0][0], ex[-1][1]))] += 1
P('| biotype | T2T-unique | SD90 | other |')
P('|---|---|---|---|')
for bt in ('protein_coding', 'lncRNA', 'other'):
    P(f'| {bt} | ' + ' | '.join(str(tab2[(bt, w)]) for w in ('T2T-unique', 'SD90', 'other')) + ' |')
P(f'CAT genes {len(CAT)}; missing in RefSeq {sum(tab2.values())}')

# ===== longer
P('\n## Same gene, CAT span >= 1.5x and >= 10 kb longer than RefSeq (coding RefSeq genes)')
maxintron = lambda trs: max((b[0] - a[1] for t in trs for a, b in zip(sorted(t), sorted(t)[1:])), default=0)
lt, rev = collections.Counter(), 0
exl = collections.defaultdict(list)
for g, (c, s, ex) in CAT.items():
    nm, src = cg[g][2], cattr.get(g, {}).get('src', '')
    for r in ri.hits(c, s, ex):
        if r not in coding:
            continue
        rn = rg[r][2]
        if not (rn == nm or (src and sym2e.get(rn) == src)):
            continue
        cs, rsp = ex[-1][1] - ex[0][0], RS[r][2][-1][1] - RS[r][2][0][0]
        if cs >= 1.5 * rsp and cs - rsp >= 10000:
            mc, mr = maxintron(ctr[g]), maxintron(rtr[r].values())
            why = ('merged' if cattr.get(g, {}).get('coll', 'N/A') not in ('', 'N/A')
                   else 'stretched intron' if mc >= 2 * mr and mc - mr >= 10000 else 'isoform extent')
            lt[why] += 1
            if len(exl[why]) < 8:
                exl[why].append(f'{nm}({cs // 1000}k vs {rsp // 1000}k)')
        elif rsp >= 1.5 * cs and rsp - cs >= 10000:
            rev += 1
for k in ('merged', 'stretched intron', 'isoform extent'):
    P(f'| {k} | {lt[k]} | ' + ' '.join(exl[k]) + ' |')
P(f'CAT longer: {sum(lt.values())}; RefSeq longer (reverse): {rev}')

# ===== merges vs RefSeq and our copies
P('\n## CAT collapsed genes: how many genes RefSeq and our assembly have at the same exons')
cp = {}
for l in open(copies):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'family_id':
        h = {k: i for i, k in enumerate(f)}
        continue
    exs = merge([(int(x.split('-')[0]), int(x.split('-')[1])) for x in f[h['exons']].split(',')])
    cp[f[h['family_id']] + ':' + f[h['copy_idx']]] = (f[h['chrom']], f[h['strand']], exs)
oi = Idx(cp)
rows = []
for g, v in cattr.items():
    if v['coll'] in ('', 'N/A') or g not in CAT:
        continue
    c, s, ex = CAT[g]
    merged = [x.split('.')[0] for x in v['coll'].split(',')]
    isp = any(m in par.get(v['src'], ()) for m in merged)
    rs_h = ri.hits(c, s, ex)
    names = set(v['colln'].split(',')) | {v['name']} | {e2sym.get(m, '') for m in merged}
    rs_named = [r for r in rs_h if rg[r][2] in names]
    rows.append((isp, len(merged) + 1, len(rs_h), len(rs_named), len(oi.hits(c, s, ex)), v['name'], c, ex[0][0], ex[-1][1]))
for lab, sel in (('paralog merges', [r for r in rows if r[0]]), ('all collapsed', rows)):
    n = len(sel)
    more_rs = sum(1 for r in sel if r[3] >= 2)
    more_us = sum(1 for r in sel if r[4] >= 2)
    P(f'{lab}: {n} CAT genes standing for {sum(r[1] for r in sel)} source genes; RefSeq has >= 2 of the merged genes at {more_rs} '
      f'({more_rs / max(n, 1):.0%}); our assembly has >= 2 copies at {more_us} ({more_us / max(n, 1):.0%})')
P('| CAT gene | locus | source genes merged | RefSeq genes (named among merged) | our copies |')
P('|---|---|---|---|---|')
for r in sorted([r for r in rows if r[0]], key=lambda r: -r[1])[:25]:
    P(f'| {r[5]} | {r[6]}:{r[7]}-{r[8]} | {r[1]} | {r[2]} ({r[3]}) | {r[4]} |')
open(out + '.md', 'w').write('\n'.join(L) + '\n')
