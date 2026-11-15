#!/usr/bin/env python3
"""Read support of gene-joining transcripts (2026-11-14, descriptive, register 1706): for each transcript of fusion_audit.py's
table, its crossing intron(s) = an intron whose left exon lies on one joined gene and whose right exon on another that
shares no exon with it; the BAM records carrying exactly that intron are counted as primary / secondary / supplementary.
Joins with no such intron (one exon over both genes, or an unannotated exon between) are counted and left out.
usage: join_support.py JOINING.tsv ASSEMBLY.gtf REFSEQ.gtf[.gz] BAM OUT.tsv     (needs pysam)"""
import collections, csv, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fusion_audit as FA
from cat_vs_refseq_lib import Idx, share


def crossing(ex, hits, R):
    """introns (start, end) of exon list `ex` between consecutive exons on two genes that share no exon."""
    on = [{g for g in hits if share([list(e)], R[g][2])} for e in ex]
    return [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)
            if any(not share(R[a][2], R[b][2]) for a in on[i] - on[i + 1] for b in on[i + 1] - on[i])]


def main(joining, gtf, rs, bam, out):
    import pysam
    R, rmeta = FA.load(rs, 'gene_id'); R = {g: v for g, v in R.items() if any(t.startswith(('rna-NM_', 'rna-XM_')) for t in rmeta[g][4])}; ri = Idx(R)
    ex = collections.defaultdict(list)
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[FA.A('transcript_id', f[8])].append((int(f[3]) - 1, int(f[4])))
    B = pysam.AlignmentFile(bam); cache = {}

    def support(c, s, e):
        if (c, s, e) not in cache:
            n = [0, 0, 0, 0]  # primary, primary MAPQ 0, secondary, supplementary
            for r in B.fetch(c, max(s - 1, 0), s):
                p = r.reference_start; hit = False
                for op, ln in r.cigartuples:
                    if op == 3 and p == s and p + ln == e: hit = True; break
                    if op in (0, 2, 3, 7, 8): p += ln
                if hit:
                    n[3 if r.is_supplementary else 2 if r.is_secondary else 0] += 1; n[1] += (not r.is_secondary and not r.is_supplementary and r.mapping_quality == 0)
            cache[(c, s, e)] = n
        return cache[(c, s, e)]

    rows = list(csv.DictReader(open(joining), delimiter='\t')); c = collections.Counter(); agg = collections.defaultdict(collections.Counter)
    with open(out, 'w') as fh:
        print('transcript\tclass\tjunc2\ttagged\treads\tintron\tprimary\tprimary_mapq0\tsecondary\tsupplementary\tsecondary_frac', file=fh)
        for r in rows:
            e = sorted(ex[r['transcript']]); X = crossing(e, ri.hits(r['chrom'], r['strand'], [list(x) for x in e]), R)
            if not X: c['no crossing intron'] += 1; continue
            s, t = min(X, key=lambda x: sum(support(r['chrom'], *x)[i] for i in (0, 2)))  # the weakest crossing intron carries the join
            p, p0, sec, sup = support(r['chrom'], s, t); tot = p + sec; fr = sec / tot if tot else float('nan')
            print(r['transcript'], r['class'], r['junc2'], r['tagged'], r['reads'], f"{r['chrom']}:{s}-{t}", p, p0, sec, sup, f'{fr:.3f}', sep='\t', file=fh)
            b = 'no read' if not tot else 'primary only' if not sec else 'primary majority' if fr < 0.5 else 'secondary majority' if p else 'secondary only'
            for k in (r['class'], 'ALL', 'ALL junc2' if r['junc2'] == '1' else 'ALL exon-overlap', 'F1-tagged' if r['tagged'] == '1' else 'untagged'):
                agg[k][b] += 1; agg[k]['n'] += 1; agg[k]['mapq0 majority'] += (p > 0 and 2 * p0 >= p)
    cols = ('n', 'primary only', 'primary majority', 'secondary majority', 'secondary only', 'no read', 'mapq0 majority')
    print(f'joining transcripts {len(rows)}; without a crossing intron {c["no crossing intron"]}; distinct crossing introns {len(cache)}')
    print('group', *cols, sep='\t')
    for k in ('ALL', 'ALL junc2', 'ALL exon-overlap', 'annot_rt', 'truth', 'named', 'loc', 'other', 'F1-tagged', 'untagged'): print(k, *(agg[k][x] for x in cols), sep='\t')


def _selfcheck():
    R = {'a': ('c', '+', [[0, 10]]), 'b': ('c', '+', [[20, 30]]), 'ab': ('c', '+', [[5, 25]])}
    assert crossing([(0, 8), (22, 30)], ['a', 'b'], R) == [(8, 22)] and crossing([(0, 25)], ['a', 'b'], R) == []
    assert crossing([(0, 8), (22, 30)], ['a', 'ab'], R) == [] and crossing([(0, 4), (6, 8), (22, 30)], ['a', 'b'], R) == [(8, 22)]


if __name__ == '__main__':
    _selfcheck(); main(*sys.argv[1:6])
