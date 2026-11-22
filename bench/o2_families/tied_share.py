#!/usr/bin/env python3
"""Read-level provenance on reads of known origin (O2 family benchmark; rows 1873, 1874, 1877).

    tied_share.py LOCI.gtf BAM AS_TABLE TRUTH.tsv OUT [THREADS]

TRUTH: read (or its parent: a name PARENT|k falls back on PARENT) <tab> locus of origin. Loci = gene_id of LOCI.gtf.
Reports reads placed alone on the true locus / another locus / no locus, the tied reads (minimap2 ms, exact band, as
tools/eq_classes.py), and on those the mean share two sharing rules give the true locus: EM (eq_classes.em) and EQUAL
(1 / class size). UNPOOL=1 gives every class its own no-locus member instead of one genome-wide `outside` (row 1877).
Writes OUT, OUT.loci.tsv (per locus: true and predicted counts), OUT.groups.tsv."""
import sys, os, collections, multiprocessing, statistics
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
import eq_classes as E
gtf, bam, tab, truthf, out = sys.argv[1:6]; thr = int(sys.argv[6]) if len(sys.argv) > 6 else 32
import pysam
meta, E.G_GEX = E.load_loci(gtf)
E.G_FLOOR = E.floor_fn(10); E.G_BAM, E.G_TAG, E.G_TIES = bam, 'ms', E.load_ties(tab, E.G_FLOOR, 'ms')
truth = {}
for l in open(truthf):
    f = l.rstrip('\n').split('\t'); truth[f[0]] = f[1]
tr = lambda q: truth.get(q) or truth.get(q.rsplit('|', 1)[0])
# placed reads need names too: rescan untied reads per contig for truth of `placed`
def scan2(c):
    pl, pairs, _ = E.scan(c)
    return pl, pairs
def wrongplaced(c):
    B = pysam.AlignmentFile(bam); from gtf_names import overlap_genes
    idx = {s: ([a for a, _, _ in E.G_GEX[(c, s)]], max(b - a for a, b, _ in E.G_GEX[(c, s)])) for s in '+-' if (c, s) in E.G_GEX}
    res = collections.Counter(); per = collections.Counter()
    for r in B.fetch(c):
        if r.is_unmapped or r.is_supplementary or r.is_secondary or r.query_name in E.G_TIES: continue
        s = E.tx_strand(r)
        ov = overlap_genes([(a + 1, b) for a, b in r.get_blocks()], E.G_GEX[(c, s)], *idx[s]) if s in idx else None
        g = max(ov, key=lambda k: (ov[k], k)) if ov else E.OUTSIDE
        t = tr(r.query_name); res['right' if g == t else 'outside' if g == E.OUTSIDE else 'wrong'] += 1; per[(t, g)] += 1
    return res, per
contigs = list(pysam.AlignmentFile(bam).references)
with multiprocessing.get_context('fork').Pool(thr) as pool:
    res = pool.map(scan2, contigs, chunksize=1); wp = pool.map(wrongplaced, contigs, chunksize=1)
base, sets = collections.Counter(), collections.defaultdict(set)
for pl, pairs in res:
    base.update(pl)
    for q, g, exact, *x in pairs:
        if exact: sets[q].add(g)
base.pop(E.OUTSIDE, None)
placed, classes, tied = collections.Counter(base), collections.Counter(), {}
WP = collections.Counter(); truep = collections.Counter()
for r, per in wp:
    WP.update(r)
    for (t, g), n in per.items(): truep[t] += n
for q, s in sets.items():
    if len(s) == 1:
        placed.update(s - {E.OUTSIDE}); g = next(iter(s)); t = tr(q); WP['right' if g == t else 'outside' if g == E.OUTSIDE else 'wrong'] += 1; truep[t] += 1
    else:
        c = frozenset(s)
        if os.environ.get('UNPOOL') and E.OUTSIDE in c: c = frozenset((c - {E.OUTSIDE}) | {'outside#' + ','.join(sorted(c - {E.OUTSIDE}))})
        classes[c] += 1; tied[q] = c
theta = E.em(classes, placed)
_pooled = collections.Counter()
for _c, _k in classes.items(): _pooled[frozenset(E.OUTSIDE if l.startswith('outside#') else l for l in _c)] += _k
loci_rows, cls = E.summarise(placed, _pooled, meta)
gid = {r[0]: r[10] for r in loci_rows}
pem = peq = 0.0; n = 0; notin = 0
G = collections.defaultdict(lambda: [0, 0.0, 0.0]); tt = collections.Counter(); pe = collections.Counter(); pq = collections.Counter()
for c, k in classes.items():
    z = sum(theta[l] for l in c)
    for l in c:
        pe[l] += k * (theta[l] / z if z > 0 else 1 / len(c)); pq[l] += k / len(c)
for q, c in tied.items():
    t = tr(q); n += 1; tt[t] += 1
    z = sum(theta[l] for l in c)
    a = (theta[t] / z if z > 0 else 1 / len(c)) if t in c else 0.0; b = 1 / len(c) if t in c else 0.0
    notin += t not in c; pem += a; peq += b
    g = gid.get(t, '.'); G[g][0] += 1; G[g][1] += a; G[g][2] += b
L = set(tt) | set(pe)
mis = lambda p: sum(abs(p[l] - tt[l]) for l in L) / (2 * max(1, n))
big = {g: v for g, v in G.items() if g != '.' and v[0] >= 50}
with open(out, 'w') as w:
    P = lambda *a: (print(*a), w.write(' '.join(map(str, a)) + '\n'))
    P('reads placed alone', sum(WP.values()), 'on true locus', WP['right'], 'on another locus', WP['wrong'], 'outside', WP['outside'])
    P('tied reads (test set)', n, 'true locus not in class', notin, 'classes', len(classes))
    P('PRIMARY mean share on true locus: EM', round(pem / max(1, n), 4), 'EQUAL', round(peq / max(1, n), 4))
    P('misallocation: EM', round(mis(pe), 4), 'EQUAL', round(mis(pq), 4))
    P('tie groups with >= 50 test reads', len(big), 'EM higher', sum(v[1] > v[2] + 1e-9 for v in big.values()), 'EQUAL higher', sum(v[2] > v[1] + 1e-9 for v in big.values()),
      'median EM', round(statistics.median([v[1] / v[0] for v in big.values()]), 4) if big else 'NA', 'median EQUAL', round(statistics.median([v[2] / v[0] for v in big.values()]), 4) if big else 'NA')
with open(out + '.loci.tsv', 'w') as w:
    w.write('locus\tgroup\ttrue_placed\ttrue_tied\tplaced\tmax\tem_tied\tequal_tied\n')
    ex = collections.Counter()
    for c, k in classes.items():
        for l in c: ex[l] += k
    for l in sorted(L | set(placed) | set(truep), key=str):
        w.write(f"{l}\t{gid.get(l, '.')}\t{truep[l]}\t{tt[l]}\t{placed.get(l, 0)}\t{placed.get(l, 0) + ex[l]}\t{pe[l]:.1f}\t{pq[l]:.1f}\n")
with open(out + '.groups.tsv', 'w') as w:
    w.write('group\ttest_reads\tem\tequal\n')
    for g, v in sorted(G.items(), key=lambda kv: -kv[1][0]): w.write(f'{g}\t{v[0]}\t{v[1] / v[0]:.4f}\t{v[2] / v[0]:.4f}\n')
