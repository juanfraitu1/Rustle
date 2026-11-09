#!/usr/bin/env python3
"""Row 1594 (docs/PREREG_rna_copynumber_soto_2026-11-09.md): copy_phasing extra-copy calls vs Soto S1E parCN classes.

  rnacn.py subset S1E.tsv CAT.tsv FAM.copies.tsv OUT.copies.tsv   families with a copy overlapping any S1E gene's exons
  rnacn.py score  S1E.tsv CAT.tsv CP_PREFIX [--per-gene OUT.tsv]  CP_PREFIX = copy_phasing --out (reads PREFIX.candidates.tsv,
                                                                  PREFIX.putative.tsv) next to cp.copies.tsv (real_net output)
  rnacn.py --selftest
"""
import sys, os, csv, collections, statistics


def exons(s):
    return sorted(tuple(map(int, x.split('-'))) for x in s.split(',') if x)


def overlap(a, b):
    """any >= 1 bp overlap between two sorted half-open-ish block lists (CAT start0-end; copies table start-end)."""
    i = j = 0
    while i < len(a) and j < len(b):
        if a[i][0] < b[j][1] and b[j][0] < a[i][1]:
            return True
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return False


def load_genes(s1e, cat):
    r = list(csv.reader(open(s1e), delimiter='\t')); h = r[0]
    ix = {c: i for i, c in reversed(list(enumerate(h)))}          # first occurrence wins (S1E repeats some names)
    loc = {}
    for x in csv.DictReader(open(cat), delimiter='\t'):
        loc[x['gene_id']] = (x['chrom'], exons(x['exon_blocks']))
    genes = []
    for x in r[1:]:
        gid, cls = x[ix['Gene ID']], x[ix['parCN class']]
        try:
            m = float(x[ix['Median parCN']])
        except ValueError:
            m = None
        lab = ('pos' if cls in ('Fixed', 'Nearly-Fixed') and m is not None and m >= 2.5 else
               'neg' if cls in ('Fixed', 'Nearly-Fixed') and m is not None and 1.5 <= m < 2.5 else
               'poly' if cls == 'Polymorphic' else 'low' if m is not None and m < 1.5 else 'other')
        if gid in loc:
            genes.append(dict(gene=x[ix['Gene Name']], gid=gid, fam=x[ix['Family ID']], cls=cls, parcn=m, label=lab,
                              chrom=loc[gid][0], exons=loc[gid][1]))
    return genes


def by_chrom(genes):
    d = collections.defaultdict(list)
    for g in genes:
        d[g['chrom']].append(g)
    return d


def copy_hits_genes(chrom, ex, gidx):
    return [g for g in gidx.get(chrom, ()) if g['exons'][0][0] < ex[-1][1] and ex[0][0] < g['exons'][-1][1] and overlap(ex, g['exons'])]


def subset(s1e, cat, copies, out):
    gidx = by_chrom(load_genes(s1e, cat))
    rows = list(csv.DictReader(open(copies), delimiter='\t'))
    keep = {r['family_id'] for r in rows if copy_hits_genes(r['chrom'], exons(r['exons']), gidx)}
    with open(out, 'w') as w:
        w.write('\t'.join(rows[0].keys()) + '\n')
        for r in rows:
            if r['family_id'] in keep:
                w.write('\t'.join(r.values()) + '\n')
    print(f'{len(keep)} families, {sum(r["family_id"] in keep for r in rows)} copies -> {out}')


def fisher_greater(a, b, c, d):
    """one-sided Fisher exact p for [[a, b], [c, d]] (rows pos/neg, cols hit/miss): P(X >= a)."""
    from math import comb
    n1, n2, k = a + b, c + d, a + c
    tot = comb(n1 + n2, k)
    return sum(comb(n1, x) * comb(n2, k - x) for x in range(a, min(n1, k) + 1)) / tot


def score(s1e, cat, prefix, per_gene=None):
    genes = load_genes(s1e, cat); gidx = by_chrom(genes)
    d = os.path.dirname(prefix) or '.'
    cps = {r['tid']: r for r in csv.DictReader(open(os.path.join(d, 'cp.copies.tsv')), delimiter='\t')}
    fam_n = collections.Counter(r['family_id'] for r in cps.values())
    testable, fam_of, reads_of = set(), {}, collections.defaultdict(int)
    for r in cps.values():
        for g in copy_hits_genes(r['chrom'], exons(r['exons']), gidx):
            testable.add(g['gid']); fam_of[g['gid']] = max(fam_of.get(g['gid'], 0), fam_n[r['family_id']])
            reads_of[g['gid']] += int(r['n_reads'])
    homes = [r['nearest_locus'] for r in csv.DictReader(open(prefix + '.candidates.tsv'), delimiter='\t') if r['flagged'] == '1']
    homes += [r['copy'] for r in csv.DictReader(open(prefix + '.putative.tsv'), delimiter='\t')]
    hit = set()
    for t in homes:
        if t in cps:
            for g in copy_hits_genes(cps[t]['chrom'], exons(cps[t]['exons']), gidx):
                hit.add(g['gid'])
    res = {}
    for lab in ('pos', 'neg', 'poly', 'low'):
        gl = [g for g in genes if g['label'] == lab]; t = [g for g in gl if g['gid'] in testable]
        h = [g for g in t if g['gid'] in hit]
        res[lab] = (len(gl), len(t), len(h))
        med = lambda v: statistics.median(v) if v else float('nan')
        print(f'{lab:5s} genes {len(gl):5d} testable {len(t):5d} hit {len(h):5d} rate {len(h) / len(t) if t else float("nan"):.3f}  '
              f'median family copies {med([fam_of[g["gid"]] for g in t])}  median copy reads {med([reads_of[g["gid"]] for g in t])}')
    (_, tp, hp), (_, tn, hn) = res['pos'], res['neg']
    if tp < 15:
        print(f'VERDICT: STOP (testable positives {tp} < 15)'); return
    sens, fpr = hp / tp, hn / tn if tn else float('nan'); p = fisher_greater(hp, tp - hp, hn, tn - hn)
    print(f'sensitivity {sens:.3f} FPR {fpr:.3f} diff {sens - fpr:+.3f} Fisher one-sided p {p:.2e} -> '
          f'{"PASS" if sens - fpr >= 0.15 and p < 0.01 else "FAIL"} (bar diff >= 0.15 and p < 0.01)')
    big = lambda lab: [g for g in genes if g['label'] == lab and g['gid'] in testable and fam_of[g['gid']] >= 3]
    bp, bn = big('pos'), big('neg'); hbp, hbn = sum(g['gid'] in hit for g in bp), sum(g['gid'] in hit for g in bn)
    if bp and bn:
        print(f'families >= 3 copies (descriptive): pos {hbp}/{len(bp)} neg {hbn}/{len(bn)} diff {hbp / len(bp) - hbn / len(bn):+.3f} '
              f'p {fisher_greater(hbp, len(bp) - hbp, hbn, len(bn) - hbn):.2e}')
    if per_gene:
        with open(per_gene, 'w') as w:
            w.write('gene\tgene_id\tsoto_family\tclass\tmedian_parCN\tlabel\ttestable\thit\tfamily_copies\tcopy_reads\n')
            for g in genes:
                w.write(f'{g["gene"]}\t{g["gid"]}\t{g["fam"]}\t{g["cls"]}\t{g["parcn"]}\t{g["label"]}\t{int(g["gid"] in testable)}\t'
                        f'{int(g["gid"] in hit)}\t{fam_of.get(g["gid"], 0)}\t{reads_of.get(g["gid"], 0)}\n')


def selftest():
    assert overlap([(0, 10)], [(9, 20)]) and not overlap([(0, 10)], [(10, 20)]) and overlap([(0, 5), (30, 40)], [(35, 36)])
    assert abs(fisher_greater(5, 0, 0, 5) - 1 / 252) < 1e-12 and abs(fisher_greater(0, 5, 5, 0) - 1.0) < 1e-12
    g = [dict(gid='G1', chrom='c', exons=[(100, 200), (300, 400)])]
    gx = by_chrom(g)
    assert copy_hits_genes('c', [(250, 290)], gx) == []                     # inside the intron: no exon overlap
    assert copy_hits_genes('c', [(390, 500)], gx) == g and copy_hits_genes('d', [(390, 500)], gx) == []
    print('selftest ok')


if __name__ == '__main__':
    a = sys.argv[1:]
    if a[0] == '--selftest':
        selftest()
    elif a[0] == 'subset':
        subset(*a[1:5])
    else:
        score(a[1], a[2], a[3], a[5] if len(a) > 5 and a[4] == '--per-gene' else None)
