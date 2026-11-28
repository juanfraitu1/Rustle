#!/usr/bin/env python3
"""Counts: does adding the unspliced reads inside a transcript's exons bring `reads` closer to the truth? (register 1989, addendum C)

    unspliced_counts.py GTF BAM TRUTH.gtf [THREADS]      simulated reads: the read name ends in the transcript id

True count = primary reads named after the transcript. Assembled transcripts are matched to the truth by exact intron chain.
An unspliced primary read wholly inside exons of several transcripts is shared in proportion to their `reads`."""
import sys, os, re, bisect, collections, multiprocessing, statistics
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread, importlib.util
_s = importlib.util.spec_from_file_location('us', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'unspliced_support.py')); US = importlib.util.module_from_spec(_s); _s.loader.exec_module(US)
_G = {}

def _scan(contig):
    import pysam
    ex, own = _G['exons'].get(contig, ([], {})); st = [e[0] for e in ex]; add = collections.Counter(); true = collections.Counter(); nu = 0
    for r in pysam.AlignmentFile(_G['bam']).fetch(contig):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        true[_G['name'](r.query_name)] += 1
        if 'N' in r.cigarstring or not ex: continue
        nu += 1; ts = {t: w for i in US.inside(st, ex, r.reference_start, r.reference_end) for t, w in own[ex[i]]}; tot = sum(ts.values())
        for t, w in ts.items(): add[t] += w / tot
    return add, true, nu

def main(gtf, bam, truth, thr=32):
    import pysam
    T = {(c, oneread.chain(ex)): t for t, (c, s, ex, _) in oneread.models(truth).items() if len(ex) > 1}
    M = {}; exons = collections.defaultdict(lambda: collections.defaultdict(list))
    for t, (c, s, ex, lines) in oneread.models(gtf).items():
        m = re.search(r' reads "([0-9.]+)"', lines[0])
        if len(ex) < 2 or not m: continue
        M[t] = (T.get((c, oneread.chain(ex))), float(m.group(1)))
        for e in ex: exons[c][e].append((t, float(m.group(1))))
    US._G['maxlen'] = max(b - a for v in exons.values() for a, b in v)
    _G.update(bam=bam, exons={c: (sorted(v), v) for c, v in exons.items()}, name=lambda q: q.rsplit('_', 1)[-1])
    with multiprocessing.get_context('fork').Pool(int(thr)) as pool: res = pool.map(_scan, list(pysam.AlignmentFile(bam).references), chunksize=1)
    add = sum((a for a, _, _ in res), collections.Counter()); true = sum((t for _, t, _ in res), collections.Counter()); nu = sum(n for _, _, n in res)
    rows = [(true[tt], r, r + add.get(t, 0)) for t, (tt, r) in M.items() if tt and true.get(tt)]
    if not rows: print('no transcript matched a truth id: check the read-name rule'); return
    err = lambda i: statistics.median(abs(x[i] - x[0]) / x[0] for x in rows)
    print(f'primary reads {sum(true.values())}, unspliced {nu} ({nu / sum(true.values()):.3f}), shared out to spliced transcripts {sum(add.values()):.0f}')
    print(f'matched transcripts {len(rows)}; median relative error of the count: reads {err(1):.4f}, reads + unspliced {err(2):.4f}; total counted / true: {sum(x[1] for x in rows) / sum(x[0] for x in rows):.4f} -> {sum(x[2] for x in rows) / sum(x[0] for x in rows):.4f}')
    big = [x for x in rows if x[2] - x[1] >= 1]; print(f'transcripts that gain >= 1 read {len(big)}: closer to the truth {sum(abs(x[2] - x[0]) < abs(x[1] - x[0]) for x in big)}, farther {sum(abs(x[2] - x[0]) > abs(x[1] - x[0]) for x in big)}')

if __name__ == '__main__': main(*sys.argv[1:])
