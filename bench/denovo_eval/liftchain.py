"""Liftoff-style exon chaining over an all-vs-all locus PAF -> extra family edges (PREREG DENOVOFAM-L, docs/LEDGER_O1.md).

  liftchain.py LOCI_GFF3 SORTED_PAF OUT_TSV [--stats STATS_TSV]

SORTED_PAF = the ungated all-vs-all (`minimap2 -x asm20 -c -X`, cg tags) with a leading column "a|b" (a < b, the two locus
names) and sorted on it, so one pair's records are adjacent (both PAF orientations). Template A lifts onto locus B when the
Liftoff chain (records on one strand, colinear in A and B, not overlapping in A; the path maximizing A-exon aligned bases)
covers >= 0.5 of A's exon union (Liftoff -a) at A-exon-weighted identity >= 0.5 (Liftoff -s) and >= 1 chained A-exon base lands
on a B exon (--min-exonic-bp 1). Edge when A lifts onto B or B onto A; w = max over directions of identity x coverage.
Records mapping genomic bases onto themselves (same contig, + strand, intervals equal within 5 bp, SELF_HIT_SLOP) are dropped.
"""
import argparse
import bisect
import re
import sys

COV, IDENT, SLOP = 0.5, 0.5, 5
CIG = re.compile(r'(\d+)([MIDNSHP=X])')


def merged(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def load_loci(path):
    """name -> (chrom, start0, local merged exon blocks [s, e))."""
    loci, cur = {}, None
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'gene':
            cur = f'{f[0]}:{f[3]}-{f[4]}'
            loci[cur] = [f[0], int(f[3]) - 1, []]
        elif f[2] == 'exon':
            loci[cur][2].append((int(f[3]) - 1, int(f[4])))
    for v in loci.values():
        s0 = v[1]
        v[2] = merged([(a - s0, b - s0) for a, b in v[2]])
    return loci


def cut(blocks, s, e):
    """Pieces of [s, e) inside the sorted disjoint blocks."""
    out = []
    i = max(0, bisect.bisect_right(blocks, (s, float('inf'))) - 1)
    while i < len(blocks) and blocks[i][0] < e:
        a, b = max(blocks[i][0], s), min(blocks[i][1], e)
        if a < b:
            out.append((a, b))
        i += 1
    return out


def overlap(blocks, s, e):
    return sum(b - a for a, b in cut(blocks, s, e))


def record_exon_bases(xb, yb, qs, qe, ts, strand, cg):
    """(X-exon bases aligned, Y-exon bases aligned, X-exon bases landing on Y exons) over the CIGAR's M/=/X columns."""
    xe = ye = both = 0
    qo = to = 0  # offsets consumed from the query (strand-aware) and the target start
    for n, op in CIG.findall(cg):
        n = int(n)
        if op in 'M=X':
            ta = ts + to
            qa = qs + qo if strand == '+' else qe - qo - n
            xe += overlap(xb, qa, qa + n)
            ye += overlap(yb, ta, ta + n)
            for u, v in cut(xb, qa, qa + n):
                if strand == '+':
                    tu, tv = ta + (u - qa), ta + (v - qa)
                else:
                    tu, tv = ta + (qa + n - v), ta + (qa + n - u)
                both += overlap(yb, tu, tv)
            qo += n
            to += n
        elif op == 'I':
            qo += n
        elif op in 'DN':
            to += n
    return xe, ye, both


def chain(recs):
    """Liftoff DAG over records (a_s, a_e, b_s, b_e, strand, a_exon, ident, landing): best total A-exon bases."""
    best = (0, 0.0, 0)
    for strand in '+-':
        rs = sorted((r for r in recs if r[4] == strand and r[5] > 0), key=lambda r: (r[0], r[1]))
        dp = []
        for j, r in enumerate(rs):
            cur = (r[5], r[5] * r[6], r[7])
            for i in range(j):
                p = rs[i]
                if p[1] > r[0]:
                    continue
                if (strand == '+' and p[3] > r[2]) or (strand == '-' and r[3] > p[2]):
                    continue
                c = dp[i]
                if c[0] + r[5] > cur[0]:
                    cur = (c[0] + r[5], c[1] + r[5] * r[6], c[2] + r[7])
            dp.append(cur)
            if cur[0] > best[0]:
                best = cur
    return best[0], (best[1] / best[0] if best[0] else 0.0), best[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gff3')
    ap.add_argument('paf')
    ap.add_argument('out')
    ap.add_argument('--stats')
    a = ap.parse_args()
    loci = load_loci(a.gff3)
    elen = {k: sum(e - s for s, e in v[2]) for k, v in loci.items()}
    w = open(a.out, 'w')
    st = open(a.stats, 'w') if a.stats else None
    if st:
        st.write('a\tb\tn_rec\tcov_ab\tid_ab\tland_ab\tcov_ba\tid_ba\tland_ba\tedge\n')
    n_pairs = n_edges = n_self = 0

    def flush(key, recs):
        nonlocal n_pairs, n_edges
        if not recs:
            return
        n_pairs += 1
        x, y = key.split('|')
        out = []
        for t, side in ((x, 0), (y, 1)):
            rr = [r[side] for r in recs]
            bases, ident, land = chain(rr)
            cov = bases / elen[t] if elen.get(t) else 0.0
            out.append((cov, ident, land, cov >= COV and ident >= IDENT and land >= 1))
        edge = max((c * i for c, i, _, ok in out if ok), default=0.0)
        if edge > 0:
            n_edges += 1
            w.write(f'{x}\t{y}\t{min(1.0, edge):.6f}\n')
        if st:
            st.write(f'{x}\t{y}\t{len(recs)}\t' + '\t'.join(f'{c:.4f}\t{i:.4f}\t{l}' for c, i, l, _ in out) + f'\t{edge:.6f}\n')

    key, recs = None, []
    for line in open(a.paf):
        f = line.rstrip('\n').split('\t')
        k, f = f[0], f[1:]
        if k != key:
            flush(key, recs)
            key, recs = k, []
        q, t = f[0], f[5]
        if q not in loci or t not in loci:
            continue
        qs, qe, strand, ts, te = int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8])
        cq, ct = loci[q], loci[t]
        if cq[0] == ct[0] and strand == '+' and abs(cq[1] + qs - ct[1] - ts) <= SLOP and abs(cq[1] + qe - ct[1] - te) <= SLOP:
            n_self += 1
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        if cg is None:
            sys.exit('PAF without cg:Z: tags')
        ident = int(f[9]) / int(f[10])
        qx, tx, both = record_exon_bases(cq[2], ct[2], qs, qe, ts, strand, cg)
        rq = (qs, qe, ts, te, strand, qx, ident, both)  # query as template
        rt = (ts, te, qs, qe, strand, tx, ident, both)  # target as template
        x, _ = k.split('|')
        recs.append((rq, rt) if q == x else (rt, rq))
    flush(key, recs)
    print(f'pairs {n_pairs}, edges {n_edges}, self-hit records dropped {n_self}', file=sys.stderr)


if __name__ == '__main__':
    main()
