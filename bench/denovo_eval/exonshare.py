#!/usr/bin/env python3
"""Row 1591 (docs/PREREG_exon_sharing_2026-11-08.md): does exon-level sharing separate true from false paralog pairs
better than exon coverage (the shipped --min-shared-exon-frac measure)? Edge level, on the families stage's gated PAF.

  exonshare.py [--species human|hg002] --gtf RUN.families.gtf [--sd SEDEF_BED] LOCI_PREFIX
     LOCI_PREFIX = PREFIX.fam.loci (reads .gff3 and .gated.paf)
  exonshare.py --locus-scores OUT.tsv [--threads N] LOCI_PREFIX   (row 1592: max SHARE per locus pair, every record)
  exonshare.py --selftest
Scores per PAF record (CIGAR cg:Z walked; locus sequences are forward-strand genomic spans):
  COV       aligned bases with an exon on both sides / exonic length of the member with the smaller exonic length
  SHARE     shared exons / exons of the member with fewer exons; an exon is shared when its projected bases overlap one
            exon of the other member by >= 50% of both lengths
  SHARE_SS  as SHARE, and every internal splice site of the exon lands within 5 bp of an internal splice site there
env: as pairscore.py (DENOVO_CAT_GFF, DENOVO_CAT_TSV, DENOVO_PARALOGS; hg002: DENOVO_HG002_GFF, DENOVO_HG002_XREF).
"""
import sys, os, re, collections, bisect
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SLACK = 5


def blocks_from_cigar(cg, qs, qe, ts, strand):
    """aligned blocks (q0, t0, n) in locus coords; for '-' the query runs backwards from qe (q0 = lowest q of the block,
    pairing q0 + n - 1 - i with t0 + i)."""
    out = []; q = qs if strand == '+' else qe; t = ts
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cg):
        n = int(n)
        if op in 'M=X':
            out.append((q if strand == '+' else q - n, t, n))
            q = q + n if strand == '+' else q - n; t += n
        elif op == 'I':
            q = q + n if strand == '+' else q - n
        elif op in 'DN':
            t += n
    return out


def seg_map(blocks, strand, from_q):
    """aligned segments as (src0, src1, dst_at_src0, dir): dst = dst_at_src0 + dir * (src - src0)."""
    out = []
    for q0, t0, n in blocks:
        if from_q:
            out.append((q0, q0 + n, t0, 1) if strand == '+' else (q0, q0 + n, t0 + n - 1, -1))
        else:
            out.append((t0, t0 + n, q0, 1) if strand == '+' else (t0, t0 + n, q0 + n - 1, -1))
    out.sort()
    return out


def project(segs, a, b):
    """dst intervals (half-open) that the aligned bases of src [a, b) map to."""
    res = []
    i = max(0, bisect.bisect_right(segs, (a,)) - 1)
    while i < len(segs) and segs[i][0] < b:
        s0, s1, d0, dr = segs[i]
        x, y = max(a, s0), min(b, s1)
        if x < y:
            u, v = d0 + dr * (x - s0), d0 + dr * (y - 1 - s0)
            res.append((min(u, v), max(u, v) + 1))
        i += 1
    return res


def pos(segs, p):
    i = bisect.bisect_right(segs, (p, float('inf'))) - 1
    if i >= 0 and segs[i][0] <= p < segs[i][1]:
        s0, _, d0, dr = segs[i]
        return d0 + dr * (p - s0)
    return None


def inter(ivs, c, d):
    return sum(max(0, min(y, d) - max(x, c)) for x, y in ivs)


def exon_bases_both(blocks, strand, ea, eb):
    """aligned bases with an exon on both sides (ea: query exons, eb: target exons; relative half-open)."""
    segs = seg_map(blocks, strand, True)
    return sum(inter(project(segs, a, b), c, d) for a, b in ea for c, d in eb)


def share(src, dst, segs):
    """(shared, shared_ss) exon counts of src projected onto dst (segs: src -> dst segments)."""
    ss_dst = sorted({b for _, b in dst[:-1]} | {a for a, _ in dst[1:]})
    near = lambda p: bool(ss_dst) and min(abs(p - s) for s in ss_dst) <= SLACK
    sh = shs = 0
    for k, (a, b) in enumerate(src):
        proj = project(segs, a, b)
        if not any(inter(proj, c, d) >= 0.5 * (b - a) and inter(proj, c, d) >= 0.5 * (d - c) for c, d in dst):
            continue
        sh += 1
        if len(src) == 1:
            shs += 1; continue
        ok = True
        for p, internal, edge in ((a, k > 0, 0), (b - 1, k < len(src) - 1, 1)):
            if not internal:
                continue
            cands = [x for x in (pos(segs, q) for q in range(p - SLACK, p + SLACK + 1)) if x is not None]
            # a splice site sits between bases: a start boundary at base a is site a, an end boundary at base b-1 is site b
            if not any(near(x) or near(x + 1) for x in cands):
                ok = False
        shs += ok
    return sh, shs


def record_scores(ea, eb, cg, qs, qe, ts, strand):
    bl = blocks_from_cigar(cg, qs, qe, ts, strand)
    la, lb = sum(b - a for a, b in ea), sum(b - a for a, b in eb)
    cov = exon_bases_both(bl, strand, ea, eb) / max(1, min(la, lb))
    if (len(ea), la) <= (len(eb), lb):
        sh, shs = share(ea, eb, seg_map(bl, strand, True)); n = len(ea)
    else:
        sh, shs = share(eb, ea, seg_map(bl, strand, False)); n = len(eb)
    return cov, sh / n, shs / n


def curve(scores, labels_true, n_truth):
    pts = []
    for k in range(21):
        th = k / 20
        tp = sum(1 for p, s in scores.items() if s >= th and p in labels_true)
        n = sum(1 for s in scores.values() if s >= th)
        pts.append((th, tp / max(1, n), tp / max(1, n_truth), n))
    return pts


def interp(pts, r):
    """precision at recall r on a curve (points sorted by recall), None if r outside the range."""
    s = sorted((p[2], p[1]) for p in pts)
    lo = [x for x in s if x[0] <= r]; hi = [x for x in s if x[0] >= r]
    if not lo or not hi:
        return None
    a, b = lo[-1], hi[0]
    return a[1] if b[0] == a[0] else a[1] + (b[1] - a[1]) * (r - a[0]) / (b[0] - a[0])


def main(argv):
    import famtruth as F, protedge as P
    opt = lambda k, d=None: argv[argv.index(k) + 1] if k in argv else d
    species, gtf, sd, prefix = opt('--species', 'human'), opt('--gtf'), opt('--sd'), argv[-1]
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
    H = P.sd_holders(loci, sd) if sd else {}
    best = collections.defaultdict(lambda: [0.0, 0.0, 0.0]); sdpair = collections.defaultdict(bool)
    nrec = 0
    for line in open(prefix + '.gated.paf'):
        f = line.rstrip('\n').split('\t')
        a, b = f[0], f[5]
        if a == b or a not in gene or b not in gene or gene[a] == gene[b] or G[gene[a]]['name'] == G[gene[b]]['name']:
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        if not cg:
            continue
        sc = record_scores(rel[a], rel[b], cg, int(f[2]), int(f[3]), int(f[7]), f[4])
        k = tuple(sorted((gene[a], gene[b]))); nrec += 1
        best[k] = [max(x, y) for x, y in zip(best[k], sc)]
        if H and P.sd_linked(H, a, b):
            sdpair[k] = True
    print(f'# {prefix}: loci {len(loci)}, mapped to expressed coding genes {len(gene)}, records scored {nrec}, '
          f'gene pairs {len(best)} (truth among them {sum(1 for k in best if k in truth)}), truth pairs {len(truth)}')
    names = ('COV', 'SHARE', 'SHARE_SS')
    curves = {nm: curve({k: v[i] for k, v in best.items()}, truth, len(truth)) for i, nm in enumerate(names)}
    print('score\tthreshold\tprecision\trecall\tpairs')
    for nm in names:
        for th, p, r, n in curves[nm]:
            print(f'{nm}\t{th:.2f}\t{p:.3f}\t{r:.3f}\t{n}')
    c60 = next(x for x in curves['COV'] if abs(x[0] - 0.60) < 1e-9)
    print(f'R0 = recall of COV >= 0.60: {c60[2]:.3f} at precision {c60[1]:.3f}')
    for nm in ('SHARE', 'SHARE_SS'):
        p = interp(curves[nm], c60[2])
        print(f'{nm}: precision at R0 {"NA (outside range)" if p is None else f"{p:.3f}"}  gain '
              f'{"NA" if p is None else f"{p - c60[1]:+.3f}"}')
    for lab, sub in (('SD-born pairs', {k for k in best if sdpair[k]}), ('not SD-born', {k for k in best if not sdpair[k]})):
        if not H or not sub:
            continue
        nt = sum(1 for k in sub if k in truth)
        cv = curve({k: best[k][0] for k in sub}, truth, nt); c = next(x for x in cv if abs(x[0] - 0.60) < 1e-9)
        out = []
        for i, nm in ((1, 'SHARE'), (2, 'SHARE_SS')):
            p = interp(curve({k: best[k][i] for k in sub}, truth, nt), c[2])
            out.append(f'{nm} gain {"NA" if p is None else f"{p - c[1]:+.3f}"}')
        print(f'{lab}: {len(sub)} pairs ({nt} truth); COV>=0.60 {c[1]:.3f}/{c[2]:.3f}; ' + ', '.join(out))
    import statistics as st
    for i, nm in enumerate(names):
        t = [v[i] for k, v in best.items() if k in truth]; fl = [v[i] for k, v in best.items() if k not in truth]
        print(f'{nm} median: true pairs {st.median(t) if t else float("nan"):.2f}, false pairs {st.median(fl) if fl else float("nan"):.2f}')


_REL = None


def _score_chunk(lines):
    out = {}
    for line in lines:
        f = line.rstrip('\n').split('\t')
        a, b = f[0], f[5]
        if a == b or a not in _REL or b not in _REL:
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        if not cg:
            continue
        sh = record_scores(_REL[a], _REL[b], cg, int(f[2]), int(f[3]), int(f[7]), f[4])[1]
        k = (min(a, b), max(a, b)); out[k] = max(out.get(k, 0.0), sh)
    return out


def locus_scores(prefix, out, threads):
    """max SHARE per locus pair over all gated-PAF records -> OUT (a b share)."""
    global _REL
    import protedge as P, multiprocessing as mp
    loci = P.read_loci(prefix + '.gff3')
    _REL = {n: [(a - (int(n.rsplit(':', 1)[1].split('-')[0]) - 1), b - (int(n.rsplit(':', 1)[1].split('-')[0]) - 1))
                for a, b in bl] for n, (_, _, bl) in loci.items()}
    lines = open(prefix + '.gated.paf').readlines()
    chunks = [lines[i:i + 20000] for i in range(0, len(lines), 20000)]
    best = {}
    with mp.Pool(threads) as pool:                     # fork: workers inherit _REL
        for part in pool.imap_unordered(_score_chunk, chunks):
            for k, v in part.items():
                best[k] = max(best.get(k, 0.0), v)
    with open(out, 'w') as w:
        for (a, b), v in sorted(best.items()):
            w.write(f'{a}\t{b}\t{v:.4f}\n')
    print(f'{len(lines)} records, {len(best)} locus pairs -> {out}')


def selftest():
    # identical 3-exon loci, '+' strand, full-length match
    ea = [(0, 100), (200, 300), (400, 500)]
    cov, sh, ss = record_scores(ea, ea, '500M', 0, 500, 0, '+')
    assert (round(cov, 3), sh, ss) == (1.0, 1.0, 1.0), (cov, sh, ss)
    # target exon 2 shifted by 20 bp: still >= 50% overlap (shared) but splice sites 20 bp off (not SS-shared)
    eb = [(0, 100), (220, 320), (400, 500)]
    cov, sh, ss = record_scores(ea, eb, '500M', 0, 500, 0, '+')
    assert sh == 1.0 and abs(ss - 2 / 3) < 1e-9 and abs(cov - 280 / 300) < 1e-9, (cov, sh, ss)
    # '-' strand: query reversed; a 500-bp locus whose exons mirror ea map onto ea exactly
    mirror = [(500 - b, 500 - a) for a, b in reversed(ea)]
    cov, sh, ss = record_scores(mirror, ea, '500M', 0, 500, 0, '-')
    assert (round(cov, 3), sh, ss) == (1.0, 1.0, 1.0), (cov, sh, ss)
    # an insertion in the query inside exon 2 of a 2-exon query; exons still shared, splice sites kept
    qa = [(0, 100), (200, 310)]; tb = [(0, 100), (200, 300)]
    cov, sh, ss = record_scores(qa, tb, '250M10I50M', 0, 310, 0, '+')
    assert sh == 1.0 and ss == 1.0, (cov, sh, ss)
    pts = [(0.0, 0.5, 0.9, 10), (0.5, 0.8, 0.6, 5), (1.0, 1.0, 0.2, 1)]
    assert abs(interp(pts, 0.75) - 0.65) < 1e-9 and interp(pts, 0.95) is None
    print('selftest ok')


if __name__ == '__main__':
    a = sys.argv[1:]
    if a == ['--selftest']:
        selftest()
    elif '--locus-scores' in a:
        locus_scores(a[-1], a[a.index('--locus-scores') + 1], int(a[a.index('--threads') + 1]) if '--threads' in a else 8)
    else:
        main(a)
