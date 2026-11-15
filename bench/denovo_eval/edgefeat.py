#!/usr/bin/env python3
"""Row 1674: de novo features of every gene pair with an alignment record (docs/PREREG_edge_rule_features_2026-11-13.md).

  edgefeat.py [--species ...] --gtf RUN.families.gtf --copies PREFIX.fam.copies.tsv --hits PAIRSCORE_HITS [--sedef BED] [--fasta SOFTMASKED_GENOME] --out TSV LOCI_PREFIX
  edgefeat.py --selftest
Pairs, stages and SHARE as edgemiss.py (a gene pair takes its locus pair with the furthest stage). Features use only the
loci, their PAF records and an optional SEDEF BED; truth and linked are labels, not features.
"""
import sys, os, collections, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exonshare as X, edgemiss as E

FEATS = ('share share_union identity aln_bp n_records exon_bp exon_frac span_cov_min span_cov_max nonexon_frac exons_min exons_max '
         'exon_len_ratio span_ratio same_chrom log_dist sd_held rank_max degree_min degree_max').split()
NEUTRAL = ['neutral_bp', 'neutral_id', 'neutral_frac']   # row 1679, with --fasta


def cover(iv, exons):
    """bases of the merged intervals iv, and of those inside exons"""
    m = E.X_merge(iv)
    return sum(e - s for s, e in m), sum(max(0, min(e, d) - max(s, c)) for s, e in m for c, d in exons)


BIN = 100_000


def sd_index(bed):
    """(chrom, 100 kb bin) -> SD mates [(start, end, mate chrom, mate start, mate end)] overlapping the bin"""
    ix = collections.defaultdict(list)
    def put(c, s, e, c2, s2, e2):
        for b in range(s // BIN, e // BIN + 1): ix[(c, b)].append((s, e, c2, s2, e2))
    for l in open(bed):
        f = l.split('\t')
        if len(f) > 5 and f[1].isdigit():
            put(f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])); put(f[3], int(f[4]), int(f[5]), f[0], int(f[1]), int(f[2]))
    return ix


def sd_held(ix, a, b, minov=300):
    """a, b = (chrom, start, end): one SD pair has a mate over each"""
    return any(min(e, a[2]) - max(s, a[1]) >= minov and c2 == b[0] and min(e2, b[2]) - max(s2, b[1]) >= minov
               for k in range(a[1] // BIN, a[2] // BIN + 1) for s, e, c2, s2, e2 in ix.get((a[0], k), ()))


COMP = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')


class Neutral:
    """per locus: sequence bytes (genome case kept) and a mask of positions that are non-exon and unmasked"""
    def __init__(self, fasta, pos, rel):
        import pysam, numpy
        self.fa, self.pos, self.rel, self.np, self.c = pysam.FastaFile(fasta), pos, rel, numpy, {}

    def get(self, n):
        if n not in self.c:
            np = self.np; c, s, e = self.pos[n]
            raw = np.frombuffer(self.fa.fetch(c, s - 1, e).encode(), dtype=np.uint8)
            ok = raw < 97                                  # upper case = not soft-masked
            for a, b in self.rel[n]: ok[max(0, a - 1):b] = False
            self.c[n] = (raw & 0xDF, ok)                   # upper-cased bases, usable mask
        return self.c[n]

    def count(self, q, t, blocks, strand):
        """(usable aligned positions, matches) over gap-free blocks (q0, t0, n) of exonshare.blocks_from_cigar"""
        np = self.np; (qs, qo), (ts, to) = self.get(q), self.get(t); bp = mt = 0
        for q0, t0, n in blocks:
            a, am = qs[q0:q0 + n], qo[q0:q0 + n]; b, bm = ts[t0:t0 + n], to[t0:t0 + n]
            if len(a) != n or len(b) != n: continue
            if strand == '-':
                a = np.frombuffer(a.tobytes().translate(COMP), dtype=np.uint8)[::-1]; am = am[::-1]
            m = am & bm; bp += int(m.sum()); mt += int((m & (a == b)).sum())
        return bp, mt


def selftest():
    import numpy as np
    nz = Neutral.__new__(Neutral); nz.np = np
    nz.c = {'q': (np.frombuffer(b'ACGTACGTAC', dtype=np.uint8), np.array([1, 1, 1, 1, 0, 0, 1, 1, 1, 1], bool)),
            't': (np.frombuffer(b'ACGAACGTAC', dtype=np.uint8), np.ones(10, bool)),
            'r': (np.frombuffer(b'GTACGTTCGT', dtype=np.uint8), np.ones(10, bool))}   # reverse complement of t
    assert nz.count('q', 't', [(0, 0, 10)], '+') == (8, 7) and nz.count('r', 't', [(0, 0, 10)], '-') == (10, 10)
    assert cover([(0, 10), (5, 20)], [(8, 12), (30, 40)]) == (20, 4)
    ix = collections.defaultdict(list); ix[('c', 0)].append((0, 1000, 'd', 0, 1000)); ix[('d', 0)].append((0, 1000, 'c', 0, 1000))
    assert sd_held(ix, ('c', 100, 900), ('d', 0, 500)) and not sd_held(ix, ('c', 100, 900), ('e', 0, 500))
    print('selftest ok')


def load_loci(prefix):
    """(loci, pos, rel) of PREFIX.gff3: name -> (chrom, strand, exons); name -> (chrom, start, end); name -> locus-relative exons"""
    import protedge as P
    loci = P.read_loci(prefix + '.gff3')
    pos = {n: (n.rsplit(':', 1)[0],) + tuple(int(x) for x in n.rsplit(':', 1)[1].split('-')) for n in loci}
    rel = {n: [(a - pos[n][1] + 1, b - pos[n][1] + 1) for a, b in bl] for n, (_, _, bl) in loci.items()}
    return loci, pos, rel


def pair_features(prefix, pos, rel, sd, keep, neu=None):
    """(a, b, stage, features) of every locus pair with a CIGAR record that keep(a, b) admits, in first-record order.
    No annotation: rank and degree count the records of ALL loci, kept or not."""
    acc, part = {}, collections.defaultdict(lambda: collections.defaultdict(float))  # part: locus -> partner -> pooled matches (all loci)
    for line in open(prefix + '.gated.paf'):
        f = line.rstrip('\n').split('\t')
        a, b = f[0], f[5]
        if a == b:
            continue
        part[a][b] += int(f[9]); part[b][a] += int(f[9])
        if not keep(a, b):
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        if not cg:
            continue
        aq = E.order(rel[a], rel[b]) if a < b else not E.order(rel[b], rel[a])
        k = (a, b) if aq else (b, a)
        if (k[1], k[0]) in acc and k not in acc:
            k, aq = (k[1], k[0]), not aq
        p = acc.setdefault(k, [E.Pair(), 0, 0, 0]); p[0].add(rel[k[0]], rel[k[1]], f, cg, aq); p[1] += 1
        if neu:
            bp, mt = neu.count(a, b, X.blocks_from_cigar(cg, int(f[2]), int(f[3]), int(f[7]), f[4]), f[4]); p[2] += bp; p[3] += mt
    rank = {n: {m: i for i, (m, _) in enumerate(sorted(d.items(), key=lambda x: -x[1]))} for n, d in part.items()}
    for (a, b), (p, nrec, nbp, nmt) in acc.items():
        st, sg, un, idy = p.stage(rel[a], rel[b])
        la, lb = sum(y - x for x, y in rel[a]), sum(y - x for x, y in rel[b])
        sa, sb = pos[a][2] - pos[a][1] + 1, pos[b][2] - pos[b][1] + 1
        ca, xa = cover([v[0] for v in p.iv], rel[a]); cb, xb = cover([v[1] for v in p.iv], rel[b])
        same = pos[a][0] == pos[b][0]
        dist = max(0, max(pos[a][1], pos[b][1]) - min(pos[a][2], pos[b][2])) if same else None
        yield a, b, st, dict(share=sg, share_union=un, identity=idy, aln_bp=p.bl, n_records=nrec, exon_bp=p.ee, exon_frac=p.ee / max(1, min(la, lb)),
                             span_cov_min=min(ca / sa, cb / sb), span_cov_max=max(ca / sa, cb / sb), nonexon_frac=1 - (xa + xb) / max(1, ca + cb),
                             exons_min=min(len(rel[a]), len(rel[b])), exons_max=max(len(rel[a]), len(rel[b])),
                             exon_len_ratio=min(la, lb) / max(la, lb), span_ratio=min(sa, sb) / max(sa, sb), same_chrom=int(same),
                             log_dist=math.log10(dist + 1) if same else 9.0, sd_held=int(sd_held(sd, pos[a], pos[b])) if sd else 0,
                             rank_max=max(rank[a][b], rank[b][a]), degree_min=min(len(part[a]), len(part[b])), degree_max=max(len(part[a]), len(part[b])),
                             neutral_bp=nbp, neutral_id=nmt / nbp if nbp >= 50 else 0.0, neutral_frac=nbp / max(1, p.bl))


def main(argv):
    import famtruth as F
    opt = lambda k, d=None: argv[argv.index(k) + 1] if k in argv else d
    species, gtf, copies, out, prefix = opt('--species', 'human'), opt('--gtf'), opt('--copies'), opt('--out'), argv[-1]
    G = F.genes_from_argv(argv, species); idx = F.exon_index(G)
    expressed = F.expressed_from_gtf(gtf, idx)
    truth, universe = F.truth_pairs(G, F.PARALOGS, F.levels('simiiformes'))
    truth = {tuple(sorted(p)) for p in truth if p[0] in expressed and p[1] in expressed and G[p[0]]['name'] != G[p[1]]['name']}
    loci, pos, rel = load_loci(prefix)
    gene = {}
    for n, (c, st, bl) in loci.items():
        g = F.map_copy(idx, c, st, bl)
        if g and G[g]['biotype'] == 'protein_coding' and g in universe and g in expressed:
            gene[n] = g
    fams, _, _ = F.load_catalog(copies, G, idx)
    fam_of = {g: k for k, v in fams.items() for g in v}
    hits = {tuple(sorted(l.split('\t')[:2])) for l in open(opt('--hits'))}
    linked = lambda k: k in hits or (fam_of.get(k[0]) is not None and fam_of.get(k[0]) == fam_of.get(k[1]))
    sd = sd_index(opt('--sedef')) if opt('--sedef') else None
    neu = Neutral(opt('--fasta'), pos, rel) if opt('--fasta') else None
    feats = FEATS + (NEUTRAL if neu else [])

    def keep(a, b):   # two expressed protein-coding genes of different symbols
        if a not in gene or b not in gene or gene[a] == gene[b] or G[gene[a]]['name'] == G[gene[b]]['name']:
            return False
        # gorilla protein-hit truth: same mapped human gene = undecidable
        return not set(G[gene[a]].get('ens') or ()) & set(G[gene[b]].get('ens') or ())
    best = {}
    for a, b, st, d in pair_features(prefix, pos, rel, sd, keep, neu):
        k = tuple(sorted((gene[a], gene[b])))
        if k in best and (E.STAGES.index(st), d['share_union']) <= (E.STAGES.index(best[k][0]), best[k][3]['share_union']):
            continue
        best[k] = (st, a, b, d)
    with open(out, 'w') as w:
        w.write('gene_a\tgene_b\tname_a\tname_b\ttruth\tlinked\tstage\t' + '\t'.join(feats) + '\tlocus_a\tlocus_b\n')
        for k, (st, la_, lb_, d) in sorted(best.items()):
            w.write(f'{k[0]}\t{k[1]}\t{G[k[0]]["name"]}\t{G[k[1]]["name"]}\t{int(k in truth)}\t{int(linked(k))}\t{st}\t'
                    + '\t'.join(f'{d[x]:.4g}' for x in feats) + f'\t{la_}\t{lb_}\n')
    print(f'{prefix}: {len(best)} gene pairs with a record; truth pairs {len(truth)}')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(sys.argv[1:])
