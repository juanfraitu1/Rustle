#!/usr/bin/env python3
"""Keep a locus only when the reads still need it under the EM (register 1905; O2 extension, joint provenance test).

    em_filter.py ASSEMBLY.gtf BAM AS_TABLE OUT_PREFIX [--min-em 2] [--label] [--pooled] [--threads N]

Loci, read-to-locus rule and ties as tools/eq_classes.py (minimap2 ms, exact band). Reads tied between loci are shared by
its EM; by default every class has its own no-locus member (row 1877), --pooled = the tool's single `outside`.
Writes OUT_PREFIX.em.tsv (locus, placed, max, em, kept) and OUT_PREFIX.gtf = the assembly without the transcripts of
loci whose EM count is under --min-em (loci no read reaches are dropped too), the others tagged reads_em.
--weighted T (register 1921, opt-in): the alignment scores enter the EM. Every alignment within 10 ms of the read's best
is a candidate with weight exp(-gap / T) (gap in ms units; a substitution costs about 5), and the EM shares the read by
abundance x weight. Reads with one candidate stay placed. `expressed` keeps its meaning (a read whose best score is at
this locus alone). T -> 0 is the default test exactly.
--label drops nothing (user 2026-11-22, the advisor's point: a family member that is found stays available): every
transcript is kept and tagged o2_status = "expressed" (reads placed here alone, EM >= --min-em), "shared" (no read placed
here alone but the EM still leaves >= --min-em: indistinguishable from its partners on these reads) or "homology_only"
(reads reach it, all better explained by other loci: the member exists, these reads do not show it expressed)."""
import sys, os, re, argparse, collections, multiprocessing
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
import eq_classes as E

def scan_w(contig):
    """as E.scan, but every candidate within the band is kept with its ms gap to the read's best: -> placed, [(q, g, gap)]"""
    import pysam
    placed, pairs = collections.Counter(), []
    idx = {s: ([x for x, _, _ in E.G_GEX[(contig, s)]], max(b - x for x, b, _ in E.G_GEX[(contig, s)])) for s in "+-" if (contig, s) in E.G_GEX}
    for r in pysam.AlignmentFile(E.G_BAM).fetch(contig):
        if r.is_unmapped or r.is_supplementary: continue
        best = E.G_TIES.get(r.query_name)
        if best is None:
            if r.is_secondary: continue
        elif r.is_secondary and r.get_tag(E.G_TAG) < E.G_FLOOR(best): continue
        s = E.tx_strand(r)
        ov = E.overlap_genes([(x + 1, b) for x, b in r.get_blocks()], E.G_GEX[(contig, s)], *idx[s]) if s in idx else None
        g = max(ov, key=lambda k: (ov[k], k)) if ov else E.OUTSIDE
        if best is None: placed[g] += 1
        else: pairs.append((r.query_name, g, 0 if not r.is_secondary else max(0, best - r.get_tag(E.G_TAG))))
    return placed, pairs

def em_w(classes, placed, T, rounds=500, tol=1e-6):
    """EM with alignment weights: classes {frozenset((locus, gap)): n}; share = theta x exp(-gap / T)."""
    import math
    cl = [([(l, math.exp(-g / T)) for l, g in c], n) for c, n in classes.items()]
    members = {l for c, _ in cl for l, _ in c}
    theta = {l: float(placed.get(l, 0)) for l in members}
    for c, n in cl:
        z = sum(w for _, w in c)
        for l, w in c: theta[l] += n * w / z
    for _ in range(rounds):
        new = {l: float(placed.get(l, 0)) for l in members}
        for c, n in cl:
            z = sum(theta[l] * w for l, w in c)
            if z > 0:
                for l, w in c: new[l] += n * theta[l] * w / z
            else:
                zz = sum(w for _, w in c)
                for l, w in c: new[l] += n * w / zz
        d = max(abs(new[l] - theta[l]) for l in members) if members else 0
        theta = new
        if d < tol: break
    return theta

def selftest():
    # two loci; 3 reads best at B by one substitution (gap 5 at A). T small: they stay at B; T = 1 and A 1000-fold
    # more abundant: A takes them (at 100-fold B still keeps about one read)
    cl = collections.Counter({frozenset({('B', 0), ('A', 5)}): 3}); pl = {'A': 3000}
    t0 = em_w(cl, pl, 1e-6); t1 = em_w(cl, pl, 1.0); t9 = em_w(cl, {'A': 30}, 1.0)
    assert abs(t0['B'] - 3) < 1e-6, t0
    assert t1['B'] < 0.5 and abs(t1['A'] + t1['B'] - 3003) < 1e-6, t1
    assert t9['B'] > 2.7, t9      # at a 10-fold abundance ratio a one-substitution lead holds (2.8 of 3 reads)
    # exact ties: weights equal, same answer as the counts-only EM
    c2 = collections.Counter({frozenset({('A', 0), ('B', 0)}): 10}); p2 = {'A': 9, 'B': 1}
    tw = em_w(c2, p2, 1.0); tc = E.em(collections.Counter({frozenset({'A', 'B'}): 10}), p2)
    assert all(abs(tw[k] - tc[k]) < 1e-4 for k in tc), (tw, tc)
    print('em_filter selftest ok')

def main():
    if sys.argv[1:] == ['--selftest']: return selftest()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('tab'); ap.add_argument('out')
    ap.add_argument('--min-em', type=float, default=2); ap.add_argument('--pooled', action='store_true'); ap.add_argument('--label', action='store_true'); ap.add_argument('--threads', type=int, default=8); ap.add_argument('--weighted', type=float, default=None, metavar='T')
    a = ap.parse_args()
    import pysam
    meta, E.G_GEX = E.load_loci(a.gtf)
    E.G_FLOOR = E.floor_fn(10); E.G_BAM, E.G_TAG, E.G_TIES = a.bam, 'ms', E.load_ties(a.tab, E.G_FLOOR, 'ms')
    with multiprocessing.get_context('fork').Pool(a.threads) as pool:
        res = pool.map(scan_w if a.weighted else E.scan, list(pysam.AlignmentFile(a.bam).references), chunksize=1)
    placed, sets = collections.Counter(), collections.defaultdict(set)
    gaps = collections.defaultdict(dict)
    for pl, pairs, *_ in res:
        placed.update(pl)
        for q, g, exact, *x in pairs:
            if a.weighted: gaps[q][g] = min(exact, gaps[q].get(g, exact)); exact = exact == 0
            if exact: sets[q].add(g)
    placed.pop(E.OUTSIDE, None); classes = collections.Counter()
    if a.weighted:
        solo = collections.Counter(); wcl = collections.Counter()      # solo = reads whose best score is at one locus alone
        for q, d in gaps.items():
            if len(d) == 1: placed.update(k for k in d if k != E.OUTSIDE); continue
            if len(sets[q]) == 1: solo.update(sets[q] - {E.OUTSIDE})
            out = 'outside#' + ','.join(sorted(k for k in d if k != E.OUTSIDE))
            wcl[frozenset((out if (k == E.OUTSIDE and not a.pooled) else k, g) for k, g in d.items())] += 1
        theta = em_w(wcl, placed, a.weighted); mx = collections.Counter(placed)
        for c, n in wcl.items():
            for l, _ in c: mx[l] += n
        for l, n in solo.items(): placed[l] += n                       # for the `expressed` label and the placed column only
        sets = {}
    for s in sets.values():
        if len(s) == 1: placed.update(s - {E.OUTSIDE})
        else:
            c = frozenset(s)
            if not a.pooled and E.OUTSIDE in c: c = frozenset((c - {E.OUTSIDE}) | {'outside#' + ','.join(sorted(c - {E.OUTSIDE}))})
            classes[c] += 1
    if not a.weighted:
        theta = E.em(classes, placed); mx = collections.Counter(placed)
        for c, n in classes.items():
            for l in c: mx[l] += n
    em = {l: theta.get(l, float(placed.get(l, 0))) for l in meta}
    keep = {l for l in meta if round(em[l], 2) >= a.min_em}   # at the precision the table prints: a 1e-4 leak to a far rival does not unseat a locus (row 1921, amendment)
    status = {l: ('expressed' if placed.get(l, 0) >= 1 else 'shared') if l in keep else 'homology_only' for l in meta}
    with open(a.out + '.em.tsv', 'w') as w:
        w.write('locus\tplaced\tmax\tem\tkept\tstatus\n')
        for l in sorted(meta, key=lambda l: (meta[l][0], meta[l][1], l)): w.write(f'{l}\t{placed.get(l, 0)}\t{mx[l]}\t{em[l]:.2f}\t{int(l in keep)}\t{status[l]}\n')
    n = 0
    with open(a.out + '.gtf.tmp', 'w') as w:
        for l in open(a.gtf):
            f = l.rstrip('\n').split('\t'); g = re.search(r'gene_id "([^"]+)"', f[8]) if len(f) > 8 else None
            if g and g.group(1) in meta and g.group(1) not in keep and not a.label: n += f[2] == 'transcript'; continue
            if g and g.group(1) in em and f[2] == 'transcript': l = l.rstrip('\n').rstrip(';') + f'; reads_em "{em[g.group(1)]:.1f}"; o2_status "{status[g.group(1)]}";\n'
            w.write(l)
    os.replace(a.out + '.gtf.tmp', a.out + '.gtf')
    c = collections.Counter(status.values())
    print(f"em_filter: {len(meta)} loci: expressed {c['expressed']}, shared {c['shared']}, homology_only {c['homology_only']} (under {a.min_em:g} reads by EM; {'kept and labelled' if a.label else str(n) + ' transcripts dropped'}) -> {a.out}.gtf")

if __name__ == '__main__':
    main()
