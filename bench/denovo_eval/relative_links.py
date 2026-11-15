#!/usr/bin/env python3
"""Relative links (register 1676, OPT-IN): locus pairs the family edge rule refused that the row-1676 score calls related.

  relative_links.py --sedef BED [--exclude-bed BED] [--model JSON] FAM_PREFIX OUT.tsv        relative_links.py --selftest
FAM_PREFIX = PREFIX.fam (reads .loci.gff3, .loci.gated.paf, .copies.tsv). No annotation, no numpy. Families are unchanged.

A link says "these two loci are related", NOT "members of one family". What was measured (docs/PREREG_relative_tier_2026-11-13.md):
on protein-coding gene pairs, one locus pair per gene pair, precision 0.830 on two gorilla samples (registered bar 0.85: not
passed); about 1/6 of the picks are young duplicates, 2/5 older paralogs, 1/3 one HGNC group only, 1/6 unrelated. The model was
fitted on A119b + HG002, so those samples are in-sample. Links between other loci (lncRNA, pseudogene, unannotated, further
locus pairs of one gene pair) were never scored against any truth.

The score is NOT a confidence outside the fitted range: training scores end at 7.2, satellite-like loci reach 30+. Column
in_fitted_range = 0 marks links whose score or any feature lies outside the training tables' range.

Population = the fitted one: a CIGAR record, loci not overlapping on one chromosome, not both in the same family, and not
passing every gate (stage != edge). --exclude-bed (the satellite BED assign uses, RUSTLE_EXCLUDE_BED) sets aside loci lying
>= 50% inside it: without it 3 of 4 links on A119b are satellite loci (row 1737), far outside what was fitted. `sd_held` is a model feature, so a SEDEF BED of this genome is required.
"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools', 'data', 'relative_links_model.json')


def scorer(m):
    logf = set(m['log1p']); z = list(zip(m['features'], m['mean'], m['sd'], m['weights']))
    return lambda d: m['intercept'] + sum(w * ((math.log1p(d[f]) if f in logf else d[f]) - mu) / sd for f, mu, sd, w in z)


def families(copies, gff):
    """locus name -> family ids. A copy row's gene_id is its locus's GFF3 id; a famreps piece is `<locus id>.pN`."""
    import re
    name = {}
    for l in open(gff):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'gene':
            name[re.search(r'Name=([^;\n]+)', f[8]).group(1)] = f'{f[0]}:{f[3]}-{f[4]}'
    fam = {}
    for i, l in enumerate(open(copies)):
        f = l.rstrip('\n').split('\t')
        if i == 0:
            h = {k: j for j, k in enumerate(f)}; continue
        g = re.sub(r'\.p\d+$', '', f[h['gene_id']])
        assert g in name, f'{copies}: copy {f[h["tid"]]} has gene_id {f[h["gene_id"]]}, which is no locus of {gff}'
        fam.setdefault(name[g], set()).add(f[h['family_id']])
    return fam


def excluded(bed, pos, frac=0.5):
    """loci with >= frac of their span inside the BED's intervals (rows named ct* are not counted)"""
    iv = {}
    for l in open(bed):
        f = l.split('\t')
        # `ct` rows of the CHM13 censat BED are centromeric-transition sequence, not satellite: NBPF, GOLGA and AGAP live there
        if len(f) > 2 and f[1].isdigit() and not (len(f) > 3 and f[3].startswith('ct')): iv.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    out = set()
    for n, (c, s, e) in pos.items():
        s -= 1; m = []
        for a, b in sorted(iv.get(c, ())):          # merged overlap, so stacked BED rows are not counted twice
            a, b = max(a, s), min(b, e)
            if a >= b: continue
            if m and a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
            else: m.append([a, b])
        if sum(b - a for a, b in m) >= frac * (e - s): out.add(n)
    return out


def selftest():
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.bed') as t:
        t.write('c\t0\t60\thsat\nc\t10\t50\nc\t200\t300\nd\t0\t10\tct_1\n'); t.flush()
        assert excluded(t.name, {'x': ('c', 1, 100), 'y': ('c', 1, 200), 'z': ('d', 1, 10)}) == {'x'}
    s = scorer(dict(features=['a', 'b'], log1p=['b'], mean=[1.0, 0.0], sd=[2.0, 1.0], weights=[2.0, 1.0], intercept=0.5))
    assert abs(s(dict(a=3.0, b=math.e - 1, c=9)) - 3.5) < 1e-12
    print('selftest ok')


def main(a):
    import edgefeat as EF
    opt = lambda k, d=None: a[a.index(k) + 1] if k in a else d
    prefix, out = a[-2], a[-1]
    if not opt('--sedef'):
        sys.exit('relative_links.py: --sedef BED is required (the score uses SD links between the two loci)')
    m = json.load(open(opt('--model', MODEL))); score = scorer(m)
    loci, pos, rel = EF.load_loci(prefix + '.loci')
    fam = families(prefix + '.copies.tsv', prefix + '.loci.gff3')

    skip = excluded(opt('--exclude-bed'), pos) if opt('--exclude-bed') else set()

    def keep(x, y):
        if x not in pos or y not in pos or x in skip or y in skip: return False
        p, q = pos[x], pos[y]
        if p[0] == q[0] and p[1] <= q[2] and q[1] <= p[2]: return False   # overlapping loci: pieces of one gene
        return not (fam.get(x) and fam.get(y) and fam[x] & fam[y])
    n = 0; rows = []
    for x, y, st, d in EF.pair_features(prefix + '.loci', pos, rel, EF.sd_index(opt('--sedef')), keep):
        n += 1
        if st == 'edge': continue
        s = score(d)
        if s >= m['threshold']:
            d['fitted'] = int(s <= m['train_score_max'] and all(lo <= d[f] <= hi for f, lo, hi in zip(m['features'], m['train_min'], m['train_max'])))
            rows.append((s, x, y, st, d))
    rows.sort(key=lambda r: (-r[0], r[1], r[2]))
    fs = lambda x: ','.join(sorted(fam[x])) if x in fam else 'NA'
    with open(out, 'w') as w:
        w.write('locus_a\tlocus_b\tfamily_a\tfamily_b\tscore\tstage\tshare\tidentity\taln_bp\tsd_held\tin_fitted_range\n')
        for s, x, y, st, d in rows:
            w.write(f"{x}\t{y}\t{fs(x)}\t{fs(y)}\t{s:.4f}\t{st}\t{d['share']:.4g}\t{d['identity']:.4g}\t{d['aln_bp']}\t{d['sd_held']}\t{d['fitted']}\n")
    both = sum(1 for r in rows if r[1] in fam and r[2] in fam); one = sum(1 for r in rows if (r[1] in fam) != (r[2] in fam))
    print(f'{len(rows)} relative links among {n} refused locus pairs with a record (score >= {m["threshold"]:.4f}): '
          f'{both} between two families, {one} family to non-family locus, {len(rows) - both - one} between non-family loci; '
          f'{sum(1 for r in rows if not r[4]["fitted"])} links outside the fitted range (in_fitted_range 0: the score is an extrapolation); '
          + (f'{len(skip)} loci set aside by --exclude-bed' if opt('--exclude-bed') else 'no --exclude-bed: satellite loci are NOT set aside'))


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(sys.argv[1:])
