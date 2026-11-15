#!/usr/bin/env python3
"""Step 0 for edge recall: the first stage of the family edge rule that drops each gene pair, truth and false alike.

  edgemiss.py [--species human|hg002|gorilla] --gtf RUN.families.gtf --copies PREFIX.fam.copies.tsv [--hits PAIRSCORE_HITS] [--out TSV] LOCI_PREFIX
  (give --hits: truth pairs are then 'linked' exactly when pairscore.py counts them; without it ~200 A119b hits read as missed)
  edgemiss.py --selftest
LOCI_PREFIX = PREFIX.fam.loci (reads .gff3 and .gated.paf). Gene pairs as exonshare.py (expressed coding genes of the
truth universe, different symbols). Per locus pair, records pooled as the binary pools them; a gene pair takes the
furthest stage any of its locus pairs reaches:
  no_record  no PAF record joins the two genes' loci (not attempted by the exon gate, or minimap2 found nothing)
  short      pooled block length < 300      low_id   pooled identity < 0.70      no_exon  no exon-to-exon base
  low_share  no single record has SHARE >= 0.5 (the shipped rule); `union` = SHARE over the union of the pair's
             records (an exon shared in any record counts once) -- the split-alignment case
  low_cov    aligned span on the locus with the longer exon union < 0.30 of that union
  edge       passes every gate: a miss here is clustering (MCL or SD damping)
env: as pairscore.py.
"""
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exonshare as X

STAGES = ('no_record', 'short', 'low_id', 'no_exon', 'low_share', 'low_cov', 'edge')


def shared_idx(src, dst, segs):
    """indices of src exons whose projection overlaps one dst exon by >= 50% of both lengths (exonshare.share's rule)."""
    out = set()
    for k, (a, b) in enumerate(src):
        proj = X.project(segs, a, b)
        if any(X.inter(proj, c, d) >= 0.5 * (b - a) and X.inter(proj, c, d) >= 0.5 * (d - c) for c, d in dst):
            out.add(k)
    return out


class Pair:
    """pooled records of one locus pair, a = the side whose exons are counted (fewer exons)."""
    __slots__ = ('nm', 'bl', 'ee', 'single', 'union', 'iv')

    def __init__(self):
        self.nm = self.bl = self.ee = 0; self.single = 0; self.union = set(); self.iv = []

    def add(self, ea, eb, f, cg, a_is_query):
        qs, qe, ts, te, st = int(f[2]), int(f[3]), int(f[7]), int(f[8]), f[4]
        bl = X.blocks_from_cigar(cg, qs, qe, ts, st)
        self.nm += int(f[9]); self.bl += int(f[10])
        q, t = (ea, eb) if a_is_query else (eb, ea)
        self.ee += X.exon_bases_both(bl, st, q, t)
        s = shared_idx(ea, eb, X.seg_map(bl, st, a_is_query))
        self.single = max(self.single, len(s)); self.union |= s
        self.iv.append(((qs, qe), (ts, te)) if a_is_query else ((ts, te), (qs, qe)))

    def stage(self, ea, eb):
        """(stage, single SHARE, union SHARE)"""
        n = len(ea); sg, un = self.single / n, len(self.union) / n
        la, lb = sum(y - x for x, y in ea), sum(y - x for x, y in eb)
        side = 0 if la >= lb else 1
        cov = sum(e - s for s, e in X_merge([v[side] for v in self.iv])) / max(1, max(la, lb))
        if self.bl < 300: st = 'short'
        elif self.nm / self.bl < 0.70: st = 'low_id'
        elif not self.ee: st = 'no_exon'
        elif sg < 0.5: st = 'low_share'
        elif cov < 0.30: st = 'low_cov'
        else: st = 'edge'
        return st, sg, un, self.nm / max(1, self.bl)


def X_merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out


def order(ea, eb):
    """True when ea is the counted side (exonshare.record_scores' choice)."""
    return (len(ea), sum(y - x for x, y in ea)) <= (len(eb), sum(y - x for x, y in eb))


def main(argv):
    import famtruth as F, protedge as P
    opt = lambda k, d=None: argv[argv.index(k) + 1] if k in argv else d
    species, gtf, copies, out, prefix = opt('--species', 'human'), opt('--gtf'), opt('--copies'), opt('--out'), argv[-1]
    G = F.genes_from_argv(argv, species); idx = F.exon_index(G)
    expressed = F.expressed_from_gtf(gtf, idx)
    truth, universe = F.truth_pairs(G, F.PARALOGS, F.levels('simiiformes'))
    truth = {tuple(sorted(p)) for p in truth if p[0] in expressed and p[1] in expressed and G[p[0]]['name'] != G[p[1]]['name']}
    loci = P.read_loci(prefix + '.gff3')
    off = {n: int(n.rsplit(':', 1)[1].split('-')[0]) - 1 for n in loci}
    rel = {n: [(a - off[n], b - off[n]) for a, b in bl] for n, (_, _, bl) in loci.items()}
    gene = {}
    for n, (c, st, bl) in loci.items():
        g = F.map_copy(idx, c, st, bl)
        if g and G[g]['biotype'] == 'protein_coding' and g in universe and g in expressed:
            gene[n] = g
    fams, _, _ = F.load_catalog(copies, G, idx)
    fam_of = {g: k for k, v in fams.items() for g in v}
    linked = lambda k: fam_of.get(k[0]) is not None and fam_of.get(k[0]) == fam_of.get(k[1])
    if opt('--hits'):   # pairscore's own hit list (PAIRSCORE_HITDUMP): the catalog mapping above under-counts linked pairs (row 1662)
        hits = {tuple(sorted(l.split('\t')[:2])) for l in open(opt('--hits'))}
        cat = linked; linked = lambda k: k in hits or cat(k)   # false pairs have no hit list: the catalog decides
    acc = {}
    for line in open(prefix + '.gated.paf'):
        f = line.rstrip('\n').split('\t')
        a, b = f[0], f[5]
        if a == b or a not in gene or b not in gene or gene[a] == gene[b] or G[gene[a]]['name'] == G[gene[b]]['name']:
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        if not cg:
            continue
        aq = order(rel[a], rel[b]) if a < b else not order(rel[b], rel[a])   # is the counted side the query?
        k = (a, b) if aq else (b, a)                                          # key: (counted locus, other locus)
        if (k[1], k[0]) in acc and k not in acc:                              # exact tie in (exons, length): one key
            k, aq = (k[1], k[0]), not aq
        acc.setdefault(k, Pair()).add(rel[k[0]], rel[k[1]], f, cg, aq)
    best = {}
    for (a, b), p in acc.items():
        s = p.stage(rel[a], rel[b]); k = tuple(sorted((gene[a], gene[b])))
        if k not in best or (STAGES.index(s[0]), s[2]) > (STAGES.index(best[k][0]), best[k][2]):
            best[k] = s
    for k in truth:
        best.setdefault(k, ('no_record', 0.0, 0.0, 0.0))
    znf = lambda k: all(G[g]['name'].startswith(('ZNF', 'ZSCAN', 'ZKSCAN', 'ZFP')) for g in k)
    rows = collections.Counter()
    for k, (s, sg, un, _) in best.items():
        cls = ('truth_linked' if linked(k) else 'truth_missed') if k in truth else ('false_linked' if linked(k) else 'false_apart')
        rows[(s, cls, 'all')] += 1
        if not znf(k): rows[(s, cls, 'noZNF')] += 1
        if s == 'low_share' and un >= 0.5:
            rows[('low_share:union>=0.5', cls, 'all')] += 1
            if not znf(k): rows[('low_share:union>=0.5', cls, 'noZNF')] += 1
    cl = ('truth_missed', 'truth_linked', 'false_apart', 'false_linked')
    print(f'# {prefix}: truth pairs {len(truth)}, gene pairs with records {len(acc)} locus pairs / {len(best)} gene pairs')
    for sub in ('all', 'noZNF'):
        print(f'[{sub}] stage\t' + '\t'.join(cl))
        for s in STAGES + ('low_share:union>=0.5',):
            print(f'[{sub}] {s}\t' + '\t'.join(str(rows[(s, c, sub)]) for c in cl))
    if out:
        with open(out, 'w') as w:
            w.write('gene_a\tgene_b\tname_a\tname_b\ttruth\tlinked\tstage\tshare\tshare_union\tidentity\n')
            for k, (s, sg, un, idy) in sorted(best.items()):
                w.write(f'{k[0]}\t{k[1]}\t{G[k[0]]["name"]}\t{G[k[1]]["name"]}\t{int(k in truth)}\t{int(linked(k))}\t{s}\t{sg:.3f}\t{un:.3f}\t{idy:.3f}\n')


def selftest():
    # three 100 bp exons each side, introns of different length: three records, one exon each
    ea = [(0, 100), (1000, 1100), (2000, 2100)]; eb = [(0, 100), (500, 600), (900, 1000)]
    p = Pair()
    for (qs, _), (ts, _) in zip(ea, eb):
        p.add(ea, eb, ['a', '', str(qs), str(qs + 100), '+', 'b', '', str(ts), str(ts + 100), '100', '100'], '100M', True)
    st, sg, un, _ = p.stage(ea, eb)
    assert (st, round(sg, 3), un) == ('low_share', 0.333, 1.0), (st, sg, un)
    q = Pair(); q.add(ea, eb, ['a', '', '0', '100', '+', 'b', '', '0', '100', '60', '100'], '100M', True)
    assert q.stage(ea, eb)[0] == 'short'
    q = Pair(); q.add(ea[:1], eb[:1], ['a', '', '0', '100', '+', 'b', '', '0', '100', '290', '300'], '100M', True)
    assert q.stage(ea[:1], eb[:1])[0] == 'edge'
    print('selftest ok')


if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main(sys.argv[1:])
