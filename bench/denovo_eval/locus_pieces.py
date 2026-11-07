#!/usr/bin/env python3
"""Row 1571: cut de novo loci into PIECES for the families input only (docs/PREREG_locus_pieces_2026-11-07.md).

In each gene_id: cut at the midpoint of every gap between merged exon blocks where both sides hold a transcript with
>= MIN_READS reads lying wholly on that side. Transcripts crossing a cut are dropped; each piece becomes gene_id GENE.pK.
The assembly GTFs are untouched; piece_merge.py undoes cuts the family stage gives no sequence reason for.

usage: locus_pieces.py IN.families.gtf OUT.gtf [--min-reads N]
       locus_pieces.py --selftest
"""
import collections, re, sys

MIN_READS = 2


def span(t):
    return min(s for s, _ in t[0]), max(e for _, e in t[0])


def cuts(tx, min_reads=MIN_READS):
    """tx: {tid: ([(s, e)], reads)} -> sorted cut positions."""
    blocks = []
    for s, e in sorted(x for t in tx.values() for x in t[0]):
        if blocks and s < blocks[-1][1]:
            blocks[-1][1] = max(blocks[-1][1], e)
        else:
            blocks.append([s, e])
    out = []
    for (_, a), (b, _) in zip(blocks, blocks[1:]):
        x = (a + b) // 2
        left = any(span(t)[1] <= x and t[1] >= min_reads for t in tx.values())
        right = any(span(t)[0] >= x and t[1] >= min_reads for t in tx.values())
        if left and right:
            out.append(x)
    return out


def pieces(tx, min_reads=MIN_READS):
    """-> (pieces as tid lists in genome order, dropped crossers)."""
    cs = cuts(tx, min_reads)
    if not cs:
        return [sorted(tx)], []
    bounds = [float('-inf')] + cs + [float('inf')]
    ps = [[t for t in sorted(tx) if lo <= span(tx[t])[0] and span(tx[t])[1] <= hi] for lo, hi in zip(bounds, bounds[1:])]
    ps = [p for p in ps if p]
    kept = {t for p in ps for t in p}
    return ps, sorted(set(tx) - kept)


def main(gtf_in, gtf_out, min_reads):
    loci = collections.defaultdict(dict)
    lines = collections.defaultdict(list)
    order = []
    for line in open(gtf_in):
        f = line.split('\t')
        if len(f) < 9:
            continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        if gid not in loci:
            order.append(gid)
        lines[(gid, tid)].append(line)
        t = loci[gid].setdefault(tid, [[], 0])
        if f[2] == 'transcript':
            r = re.search(r'\breads "([^"]+)"', f[8])
            t[1] = float(r.group(1)) if r else 0
        elif f[2] == 'exon':
            t[0].append((int(f[3]) - 1, int(f[4])))
    n_split = n_pieces = n_drop = 0
    with open(gtf_out, 'w') as oh:
        for gid in order:
            tx = loci[gid]
            ps, dropped = pieces(tx, min_reads) if len(tx) >= 2 else ([sorted(tx)], [])
            if len(ps) == 1 and not dropped:
                for t in ps[0]:
                    oh.writelines(lines[(gid, t)])
                continue
            n_split += 1
            n_pieces += len(ps)
            n_drop += len(dropped)
            for k, p in enumerate(ps, 1):
                for t in p:
                    oh.writelines(l.replace(f'gene_id "{gid}"', f'gene_id "{gid}.p{k}"') for l in lines[(gid, t)])
    print(f"[locus_pieces] loci {len(loci)}, split {n_split} into {n_pieces} pieces, crossers dropped {n_drop}", file=sys.stderr)


def selftest():
    # two genes A (0-300) and B (1000-1300) joined by readthroughs R1, R2: one cut, R1/R2 dropped
    tx = {'A1': ([(0, 100), (200, 300)], 5), 'B1': ([(1000, 1100), (1200, 1300)], 4),
          'R1': ([(0, 100), (200, 300), (1000, 1100)], 3), 'R2': ([(200, 300), (1000, 1100)], 9)}
    ps, dr = pieces(tx)
    assert ps == [['A1'], ['B1']] and dr == ['R1', 'R2'], (ps, dr)
    # B supported by a 1-read transcript only: no cut
    tx['B1'] = ([(1000, 1100), (1200, 1300)], 1)
    assert pieces(tx) == ([sorted(tx)], [])
    # isoforms of ONE gene in a chain are cut too (A2 crosses, A1/A3 each lie on a side): read structure cannot tell this
    # from two genes (step 0); piece_merge.py rejoins pieces that land in one family without aligning to each other
    tx = {'A1': ([(0, 100), (200, 300)], 5), 'A2': ([(200, 300), (400, 500)], 5), 'A3': ([(400, 500), (600, 700)], 5)}
    assert pieces(tx) == ([['A1'], ['A3']], ['A2'])
    # a gap no side transcript clears (a 2-exon gene's intron) is never cut
    assert pieces({'X': ([(0, 100), (500, 600)], 9)}) == ([['X']], [])
    print('selftest ok')


if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        selftest()
    else:
        a = sys.argv[1:]
        main(a[0], a[1], float(a[a.index('--min-reads') + 1]) if '--min-reads' in a else MIN_READS)
