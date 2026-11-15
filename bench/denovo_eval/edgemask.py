#!/usr/bin/env python3
"""PREREG row 1693, step 0: lowercase share of the exon-to-exon aligned bases under each `edge`-stage gene pair.

  edgemask.py FASTA [edgemiss.py arguments ...] LOCI_PREFIX     (arguments as edgemiss.py; --out TSV writes the pairs)
  edgemask.py --selftest"""
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exonshare as X, edgemiss as E


def both(blocks, strand, eq, et):
    """(query intervals, target intervals) of the aligned bases that are exonic on both sides, locus coords, half-open."""
    Q = []; T = []
    for q0, t0, n in blocks:
        for a, b in eq:
            x, y = max(a, q0), min(b, q0 + n)
            if x >= y: continue
            u0, v0 = (t0 + x - q0, t0 + y - q0) if strand == '+' else (t0 + q0 + n - y, t0 + q0 + n - x)
            for c, d in et:
                u, v = max(c, u0), min(d, v0)
                if u >= v: continue
                T.append((u, v)); Q.append((u - t0 + q0, v - t0 + q0) if strand == '+' else (q0 + n - (v - t0), q0 + n - (u - t0)))
    return Q, T


def low_share(fa, spans):
    """spans: (chrom, start, end) 0-based half-open -> (lowercase, total)."""
    lo = tot = 0
    for c, s, e in spans:
        q = fa.fetch(c, s, e); tot += len(q); lo += sum(1 for ch in q if ch.islower())
    return lo, tot


def main(argv):
    import pysam, famtruth as F, protedge as P
    fa = pysam.FastaFile(argv[0]); argv = argv[1:]
    opt = lambda k, d=None: argv[argv.index(k) + 1] if k in argv else d
    species, gtf, copies, out, prefix = opt('--species', 'human'), opt('--gtf'), opt('--copies'), opt('--out'), argv[-1]
    G = F.genes_from_argv(argv, species); idx = F.exon_index(G); expressed = F.expressed_from_gtf(gtf, idx)
    truth, universe = F.truth_pairs(G, F.PARALOGS, F.levels('simiiformes'))
    truth = {tuple(sorted(p)) for p in truth if p[0] in expressed and p[1] in expressed and G[p[0]]['name'] != G[p[1]]['name']}
    loci = P.read_loci(prefix + '.gff3')
    off = {n: int(n.rsplit(':', 1)[1].split('-')[0]) - 1 for n in loci}; chrom = {n: n.rsplit(':', 1)[0] for n in loci}
    rel = {n: [(a - off[n], b - off[n]) for a, b in bl] for n, (_, _, bl) in loci.items()}
    gene = {}
    for n, (c, st, bl) in loci.items():
        g = F.map_copy(idx, c, st, bl)
        if g and G[g]['biotype'] == 'protein_coding' and g in universe and g in expressed: gene[n] = g
    fams, _, _ = F.load_catalog(copies, G, idx); fam_of = {g: k for k, v in fams.items() for g in v}
    hits = {tuple(sorted(l.split('\t')[:2])) for l in open(opt('--hits'))} if opt('--hits') else set()
    linked = lambda k: k in hits or (fam_of.get(k[0]) is not None and fam_of.get(k[0]) == fam_of.get(k[1]))
    acc = {}; spans = collections.defaultdict(list)
    for line in open(prefix + '.gated.paf'):   # the pooling below is edgemiss.main's, plus the genome spans of the shared exon bases
        f = line.rstrip('\n').split('\t'); a, b = f[0], f[5]
        if a == b or a not in gene or b not in gene or gene[a] == gene[b] or G[gene[a]]['name'] == G[gene[b]]['name']: continue
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        if not cg: continue
        aq = E.order(rel[a], rel[b]) if a < b else not E.order(rel[b], rel[a]); k = (a, b) if aq else (b, a)
        if (k[1], k[0]) in acc and k not in acc: k, aq = (k[1], k[0]), not aq
        acc.setdefault(k, E.Pair()).add(rel[k[0]], rel[k[1]], f, cg, aq)
        Q, T = both(X.blocks_from_cigar(cg, int(f[2]), int(f[3]), int(f[7]), f[4]), f[4], rel[a], rel[b])
        spans[k] += [(chrom[a], off[a] + x, off[a] + y) for x, y in Q] + [(chrom[b], off[b] + x, off[b] + y) for x, y in T]
    best = {}
    for (a, b), p in acc.items():
        s = p.stage(rel[a], rel[b]); k = tuple(sorted((gene[a], gene[b])))
        if k not in best or (E.STAGES.index(s[0]), s[2]) > (E.STAGES.index(best[k][0][0]), best[k][0][2]): best[k] = (s, (a, b))
    rows = []
    for k, (s, lp) in best.items():
        if s[0] != 'edge' or not linked(k): continue
        lo, tot = low_share(fa, spans[lp]); rows.append((k, int(k in truth), lo / max(1, tot), tot))
    for t, name in ((0, 'false linked'), (1, 'true linked')):
        v = [r[2] for r in rows if r[1] == t]; n = max(1, len(v))
        print(f'{name} edge-stage pairs {len(v)}: LOW >= 0.9 {sum(x >= 0.9 for x in v)} ({sum(x >= 0.9 for x in v) / n:.3f}); LOW >= 0.5 {sum(x >= 0.5 for x in v)} ({sum(x >= 0.5 for x in v) / n:.3f}); '
              + 'deciles ' + ' '.join(f'{sorted(v)[min(len(v) - 1, int(q * len(v) / 10))]:.2f}' for q in range(0, 11)) if v else f'{name}: none')
    if out:
        with open(out, 'w') as w:
            w.write('gene_a\tgene_b\tname_a\tname_b\ttruth\tlow\tbases\n')
            for k, t, lo, tot in sorted(rows): w.write(f'{k[0]}\t{k[1]}\t{G[k[0]]["name"]}\t{G[k[1]]["name"]}\t{t}\t{lo:.3f}\t{tot}\n')


def selftest():
    # plus strand: one 100-base block, query exon 10-60, target exon 30-90 (shifted by +5): shared = q 25-60 / t 30-65
    assert both([(0, 5, 100)], '+', [(10, 60)], [(30, 90)]) == ([(25, 60)], [(30, 65)])
    # minus strand: q0 + n - 1 - i pairs with t0 + i; q 0-10 pairs with t 90-100, of which the target exon keeps 95-100 = q 0-5
    assert both([(0, 0, 100)], '-', [(0, 10)], [(95, 200)]) == ([(0, 5)], [(95, 100)])
    assert both([(0, 0, 100)], '+', [(0, 10)], [(50, 60)]) == ([], [])
    print('selftest ok')


if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main(sys.argv[1:])
