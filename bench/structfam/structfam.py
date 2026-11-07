#!/usr/bin/env python3
"""Row 1569: sequence-free structure view of O1 families (docs/PREREG_structure_view_2026-11-07.md).

usage: structfam.py RUN_PREFIX OUT_TSV      (RUN_PREFIX e.g. ~/runs/human_A119b/a1)
       structfam.py --selftest
"""
import bisect, collections, itertools, math, random, sys

TOL, ALIKE, NULL_K, CAP, SEED = 3, 0.5, 20, 12, 1569


def lcs_match(a, b, tol=TOL):
    """Longest in-order chain of exon lengths matching within tol bp (skips allowed on both sides)."""
    prev = [0] * (len(b) + 1)
    for x in a:
        cur = [0]
        for j, y in enumerate(b):
            cur.append(prev[j] + 1 if abs(x - y) <= tol else max(prev[j + 1], cur[j]))
        prev = cur
    return prev[-1]


def alike(a, b):
    return lcs_match(a, b) / min(len(a), len(b)) >= ALIKE


def internal(exons, strand):
    """exons: [(start, end_exclusive)] -> internal exon lengths 5'->3' (None when < 3 exons)."""
    ex = sorted(exons)
    if len(ex) < 3:
        return None
    lens = [e - s for s, e in ex]
    if strand == '-':
        lens.reverse()
    return tuple(lens[1:-1])


def overlap(a, b):
    return a['chrom'] == b['chrom'] and a['start'] < b['end'] and b['start'] < a['end']


def fisher_greater(a, b, c, d):
    """One-sided p that row 1 (a alike / b not) is enriched over row 2 (c / d): hypergeometric upper tail."""
    n1, k, n = a + b, a + c, a + b + c + d
    lf = lambda x: math.lgamma(x + 1)
    def pmf(x):
        return math.exp(lf(k) - lf(x) - lf(k - x) + lf(n - k) - lf(n1 - x) - lf(n - k - n1 + x) - lf(n) + lf(n1) + lf(n - n1))
    return min(1.0, sum(pmf(x) for x in range(a, min(n1, k) + 1)))


def load(prefix):
    members = []
    with open(prefix + '.fam.copies.tsv') as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            r = dict(zip(hdr, line.rstrip('\n').split('\t')))
            ex = [tuple(map(int, x.split('-'))) for x in r['exons'].split(',')]
            members.append(dict(fam=r['family_id'], tid=r['tid'], chrom=r['chrom'], start=int(r['locus_start']),
                                end=int(r['locus_end']), n_exon=len(ex), ints=internal(ex, r['strand']),
                                key=f"{r['chrom']}:{int(r['locus_start']) + 1}-{r['locus_end']}",
                                gene_id=r['gene_id'], strand=r['strand'], exons=r['exons']))
    spans = collections.defaultdict(list)
    for m in members:
        spans[m['chrom']].append((m['start'], m['end']))
    for v in spans.values():
        v.sort()
    maxlen = {c: max(e - s for s, e in v) for c, v in spans.items()}

    def hits_member(chrom, s, e):
        v = spans.get(chrom)
        if not v:
            return False
        i = bisect.bisect_left(v, (s - maxlen[chrom], -1))
        return any(ms < e and s < me for ms, me in v[i:bisect.bisect_left(v, (e, -1))])

    genes, exons = {}, collections.defaultdict(list)
    with open(prefix + '.fam.loci.gff3') as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) < 9:
                continue
            if f[2] == 'gene':
                gid = f[8].split('ID=gene-')[1].split(';')[0].strip()
                genes[gid] = dict(chrom=f[0], start=int(f[3]) - 1, end=int(f[4]), strand=f[6])
            elif f[2] == 'exon':
                exons[f[8].split('gene=')[1].strip()].append((int(f[3]) - 1, int(f[4])))
    null = []
    for gid, g in genes.items():
        if hits_member(g['chrom'], g['start'], g['end']):
            continue
        ints = internal(exons[gid], g['strand'])
        if ints:
            null.append(dict(g, fam=None, tid=gid, gene_id=gid, ints=ints,
                             exons=','.join(f'{s}-{e}' for s, e in sorted(exons[gid]))))
    return members, null, len(genes)


def pair_identity(prefix, keys):
    ident = {}
    with open(prefix + '.fam.loci.gated.paf') as fh:
        for line in fh:
            f = line.split('\t', 12)
            if f[0] in keys and f[5] in keys and f[0] != f[5]:
                k = frozenset((f[0], f[5]))
                ident[k] = max(ident.get(k, 0.0), int(f[9]) / int(f[10]))
    return ident


def matched_null(pairs, buckets, rng):
    """NULL_K n_int-matched non-overlapping pairs per pair; returns (alike, total, skipped pairs)."""
    hit = tot = skipped = 0
    for a, b in pairs:
        ba, bb = buckets.get(min(len(a['ints']), CAP)), buckets.get(min(len(b['ints']), CAP))
        if not ba or not bb:
            skipped += 1
            continue
        for _ in range(NULL_K):
            for _try in range(50):
                x, y = rng.choice(ba), rng.choice(bb)
                if x is not y and not overlap(x, y):
                    hit += alike(x['ints'], y['ints'])
                    tot += 1
                    break
    return hit, tot, skipped


def bar(name, pairs, buckets, rng, out):
    a = sum(alike(x['ints'], y['ints']) for x, y in pairs)
    nh, nt, sk = matched_null(pairs, buckets, rng)
    rf, rn = a / max(len(pairs), 1), nh / max(nt, 1)
    p = fisher_greater(a, len(pairs) - a, nh, nt - nh)
    ok = a > 0 and (rn == 0 or rf >= 3 * rn) and p < 1e-6
    out.append((name, f"pairs={len(pairs)} alike={a} rate={rf:.4f} | null pairs={nt} alike={nh} rate={rn:.4f} "
                      f"(skipped {sk}) | ratio={rf / rn if rn else float('inf'):.2f} p={p:.3g}", 'PASS' if ok else 'FAIL'))
    return ok


def main(prefix, out_tsv):
    rng = random.Random(SEED)
    members, null, n_loci = load(prefix)
    test = [m for m in members if m['ints']]
    fams = collections.defaultdict(list)
    for m in test:
        fams[m['fam']].append(m)
    F = [(a, b) for v in fams.values() for a, b in itertools.combinations(v, 2) if not overlap(a, b)]
    ident = pair_identity(prefix, {m['key'] for m in test})
    idf = lambda a, b: ident.get(frozenset((a['key'], b['key'])))
    F_lo = [(a, b) for a, b in F if (i := idf(a, b)) is not None and i < 0.95]
    buckets = collections.defaultdict(list)
    for x in null:
        buckets[min(len(x['ints']), CAP)].append(x)

    out = [('substrate', f"loci={n_loci} members={len(members)} testable_members={len(test)} "
                         f"null_testable={len(null)} families_with_2plus_testable={sum(len(v) >= 2 for v in fams.values())}", '')]
    b1 = bar('B1 F', F, buckets, rng, out)
    b2 = bar('B2 F_lo', F_lo, buckets, rng, out)

    # B3: exact internal-exon tuples with n_int >= 3, members + null universe
    groups = collections.defaultdict(list)
    for x in test + null:
        if len(x['ints']) >= 3:
            groups[x['ints']].append(x)
    same = tot = 0
    for v in groups.values():
        for a, b in itertools.combinations(v, 2):
            if not overlap(a, b):
                tot += 1
                same += a['fam'] is not None and a['fam'] == b['fam']
    prec = same / max(tot, 1)
    o1 = [(a, b) for a, b in F if len(a['ints']) >= 3 and len(b['ints']) >= 3]
    rec = sum(a['ints'] == b['ints'] for a, b in o1) / max(len(o1), 1)
    out.append(('B3 exact tuples', f"pairs={tot} same_O1_family={same} precision={prec:.4f}", 'PASS' if prec >= 0.5 else 'FAIL'))
    out.append(('desc B3 recall', f"O1 pairs n_int>=3: {len(o1)}, identical tuple {rec:.4f}", ''))

    # descriptive: identity bins
    bins = [(0, .80), (.80, .90), (.90, .95), (.95, .99), (.99, 1.01)]
    for lo, hi in bins:
        sub = [(a, b) for a, b in F if (i := idf(a, b)) is not None and lo <= i < hi]
        k = sum(alike(a['ints'], b['ints']) for a, b in sub)
        out.append((f'desc identity [{lo:.2f},{hi:.2f})', f"pairs={len(sub)} alike={k} rate={k / max(len(sub), 1):.4f}", ''))
    nopaf = [(a, b) for a, b in F if idf(a, b) is None]
    k = sum(alike(a['ints'], b['ints']) for a, b in nopaf)
    out.append(('desc identity none (MCL-transitive)', f"pairs={len(nopaf)} alike={k} rate={k / max(len(nopaf), 1):.4f}", ''))

    # descriptive: "all members alike" families vs matched null families
    big = [v for v in fams.values() if len(v) >= 3]
    def all_alike(v):
        return all(alike(a['ints'], b['ints']) for a, b in itertools.combinations(v, 2) if not overlap(a, b))
    obs = sum(all_alike(v) for v in big)
    nul = 0
    for v in big:
        for _ in range(NULL_K):
            nul += all_alike([rng.choice(buckets[min(len(m['ints']), CAP)]) for m in v])
    out.append(('desc all-alike families', f"families>=3 testable={len(big)} all_alike={obs} rate={obs / max(len(big), 1):.4f} "
                                         f"| null rate={nul / max(NULL_K * len(big), 1):.4f}", ''))
    verdict = 'PASS' if b1 and b2 else 'FAIL'
    out.append(('VERDICT B1+B2', verdict, ''))
    with open(out_tsv, 'w') as fh:
        for r in out:
            fh.write('\t'.join(r) + '\n')
    for r in out:
        print('\t'.join(r))


def selftest():
    assert lcs_match((100, 200, 300), (100, 203, 300)) == 3
    assert lcs_match((100, 200, 300), (100, 204, 300)) == 2
    assert lcs_match((100, 300), (100, 200, 300)) == 2            # skip on one side
    assert lcs_match((300, 100), (100, 300)) == 1                 # order matters
    assert alike((100, 200), (100, 999)) and not alike((100, 200, 300), (999, 998, 300))
    assert internal([(0, 10), (20, 50), (60, 100)], '+') == (30,)
    assert internal([(0, 10), (20, 50), (60, 100), (200, 260)], '-') == (40, 30)
    assert internal([(0, 10), (20, 50)], '+') is None
    assert abs(fisher_greater(10, 0, 0, 10) - 1 / 184756) < 1e-12
    print('selftest ok')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:3])
