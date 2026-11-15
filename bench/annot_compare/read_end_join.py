#!/usr/bin/env python3
"""Row 1713 (docs/PREREG_read_end_join_tag_2026-11-14.md): gene-join tag from read ends. A transcript is tagged when one
of its introns holds a site where reads coming from upstream end, followed by a site where reads running downstream start.
Scored by blind_join.main (row 1707's bars). ENDS = directory of per-chromosome start / end arrays of primary records.
usage: read_end_join.py ends BAM ENDS_DIR THREADS
       read_end_join.py X.gtf REFSEQ.gtf[.gz] JOINS.tsv ENDS_DIR GENOME.fa OUTPREFIX [TANDEM_FUSED.tsv]
writes OUTPREFIX.<form>.list.tsv (tid, chrom, strand, junctions, label, gene: a RUSTLE_BRIDGE_UNITS_LIST)."""
import collections, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blind_join

MARGIN, WIN, ARICH = 10, 50, 12


def table(chrom, strand, start, end):
    """per (chrom, strand): 5' and 3' positions in transcript orientation (minus strand negated), sorted both ways."""
    D = {}
    for c in np.unique(chrom):
        for st in '+-':
            m = (chrom == c) & (strand == st)
            if not m.any(): continue
            f, t = (start[m], end[m]) if st == '+' else (-end[m], -start[m]); i, j = np.argsort(t, kind='stable'), np.argsort(f, kind='stable')
            D[(c, st)] = (t[i], f[i], f[j], t[j])
    return D


def sites(x, n):
    """first positions of WIN-bp windows holding >= n of the sorted positions x."""
    return x[:len(x) - n + 1][x[n - 1:] - x[:len(x) - n + 1] < WIN] if len(x) >= n else x[:0]


def ends_starts(D, key, d, a):
    if key not in D or a - d <= 2 * MARGIN: return np.empty(0, int), np.empty(0, int)
    t3, f5, s5, e3 = D[key]
    i, j = np.searchsorted(t3, d + MARGIN), np.searchsorted(t3, a - MARGIN, 'right'); p = t3[i:j][f5[i:j] < d]
    i, j = np.searchsorted(s5, d + MARGIN), np.searchsorted(s5, a - MARGIN, 'right'); q = s5[i:j][e3[i:j] > a]
    return p, np.sort(q)


def arich(fa, key, p):
    if key[1] == '+': return fa.fetch(key[0], p, p + 20).upper().count('A') >= ARICH
    return fa.fetch(key[0], max(-p - 20, 0), -p).upper().count('T') >= ARICH


def signature(p, q, n, fa=None, key=None):
    """a site of >= n end reads upstream of a site of >= n start reads."""
    sq = sites(q, n)
    if not len(sq): return False
    for x in sites(p, n):
        if x >= sq[-1]: return False
        if fa is None or not arich(fa, key, int(x)): return True
    return False


def tagger(D, fa, out, pre=''):
    def run(X):
        M = table(*(np.array(v) for v in zip(*((v['c'], v['s'], v['ex'][0][0], v['ex'][-1][1]) for v in X.values()))))
        tags = collections.defaultdict(set); lists = collections.defaultdict(dict); cache = {}; t0 = time.time()
        jr = collections.Counter()  # reads of all transcripts carrying an intron (form E5)
        for v in X.values():
            for i in v['chain'][2]: jr[(v['c'], v['s'], i)] += v['reads']
        for k, (t, v) in enumerate(X.items()):
            if k % 20000 == 0: print(f'[read_end] {k} transcripts {time.time() - t0:.0f} s', file=sys.stderr, flush=True)
            key = (v['c'], v['s']); r = max(5, v['reads'])
            for s, e in v['chain'][2]:
                d, a = (s, e) if v['s'] == '+' else (-e, -s)
                if (key, d, a) not in cache:
                    p, q = ends_starts(D, key, d, a); mp, mq = ends_starts(M, key, d, a)
                    e2 = signature(p, q, 5, fa, key)
                    cache[(key, d, a)] = (p, q, {'ceiling N=1': signature(p, q, 1), 'E1': signature(p, q, 5), 'E2': e2, 'E4': e2 and signature(p, q, 10, fa, key), 'M models': signature(mp, mq, 1)})
                p, q, f = cache[(key, d, a)]
                r5 = max(5, jr[(v['c'], v['s'], (s, e))])
                for n, y in list(f.items()) + [('E3', f['E2'] and (r == 5 or signature(p, q, r, fa, key))), ('E5', f['E2'] and (r5 == 5 or signature(p, q, r5, fa, key)))]:
                    if y: tags[n].add(t); lists[n].setdefault(t, []).append(f'{s + 1}-{e}')  # 1-based intron start, as x.bridge_junctions.tsv
        for n, L in lists.items():
            if n[0] != 'E': continue
            with open(f'{out}.{n}.list.tsv', 'w') as fh:
                print('tid\tchrom\tstrand\tjunctions\tlabel\tgene', file=fh)
                for t, js in L.items(): print(t, X[t]['c'], X[t]['s'], ';'.join(js), 'readend', X[t]['g'], sep='\t', file=fh)
        return {n: tags[n] for n in ('ceiling N=1', 'M models', 'E1', 'E2', 'E3', 'E4', 'E5')}
    return run


def _ends_one(job):
    import pysam
    bam, c, out = job; st, a, b = [], [], []
    for r in pysam.AlignmentFile(bam).fetch(c):
        if r.flag & 2308: continue  # primary records only, no MAPQ filter
        rev = r.is_reverse
        if r.has_tag('ts') and r.get_tag('ts') == '-': rev = not rev  # transcript strand
        st.append(rev); a.append(r.reference_start); b.append(r.reference_end)
    np.savez(os.path.join(out, c + '.npz'), minus=np.array(st, bool), start=np.array(a, np.int64), end=np.array(b, np.int64)); return c, len(a)


def ends(bam, out, threads):
    import multiprocessing, pysam
    os.makedirs(out, exist_ok=True); cs = pysam.AlignmentFile(bam).references
    with multiprocessing.Pool(int(threads)) as p:
        for c, n in p.imap_unordered(_ends_one, [(bam, c, out) for c in cs]): print(c, n, flush=True)


def main(xgtf, rs, joins, ends_dir, fasta, out, fused=None):
    import pysam
    D = {}
    for f in sorted(os.listdir(ends_dir)):
        z = np.load(os.path.join(ends_dir, f)); c = f[:-4]
        D.update(table(np.full(len(z['start']), c), np.where(z['minus'], '-', '+'), z['start'], z['end']))
    print(f'[read_end] {sum(len(v[0]) for v in D.values())} primary records', file=sys.stderr, flush=True)
    blind_join.main(xgtf, rs, joins, fused, tagger(D, pysam.FastaFile(fasta), out))


def _selfcheck():
    assert list(sites(np.array([1, 2, 3, 100, 400]), 3)) == [1] and not len(sites(np.array([1, 100]), 2))
    # plus strand, intron (100, 1000): 5 reads from upstream end near 300, 5 reads start near 600 and run past 1000
    c = np.array(['c'] * 10); st = np.array(['+'] * 10); D = table(c, st, np.array([0] * 5 + [600, 601, 602, 603, 604]), np.array([300, 301, 302, 303, 304] + [2000] * 5))
    p, q = ends_starts(D, ('c', '+'), 100, 1000); assert len(p) == 5 and len(q) == 5 and signature(p, q, 5) and not signature(p, q, 6)
    D = table(c, st, np.array([0] * 5 + [200] * 5), np.array([700] * 5 + [2000] * 5)); assert not signature(*ends_starts(D, ('c', '+'), 100, 1000), 5)  # start before end: one gene's isoforms
    D = table(c, np.array(['-'] * 10), np.array([700] * 5 + [0] * 5), np.array([2000] * 5 + [400] * 5)); assert signature(*ends_starts(D, ('c', '-'), -1000, -100), 5)  # mirrored


if __name__ == '__main__':
    _selfcheck(); ends(*sys.argv[2:5]) if sys.argv[1] == 'ends' else main(*sys.argv[1:8])
