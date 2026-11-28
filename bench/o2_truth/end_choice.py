#!/usr/bin/env python3
"""Which read end should a transcript take? Candidates from the reads of each chain, against several judges (register 1992).

    end_choice.py GTF BAM TRUTH.gtf|- OUT.tsv [--cage PEAKS.bed] [--polya SITES.bed] [--threads N]

Chains = the multi-exon transcripts of GTF (without one-read tags). Reads = primary alignments with exactly the chain.
Prints the table; OUT.tsv = one line per chain with every candidate position."""
import sys, os, bisect, collections, multiprocessing, statistics
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread, one_read_loci as OL
W = 10
_G = {}

def candidates(v, up):
    """v = read end positions; up = +1 when a larger coordinate is the more extended end, -1 otherwise -> dict of positions"""
    s = sorted(v, reverse=(up > 0)); n = len(s); out = {'ext': s[0], 'med': sorted(v)[(n - 1) // 2] if up < 0 else sorted(v)[n // 2]}
    a = sorted(v); best = max(range(n), key=lambda i: (bisect.bisect_right(a, a[i] + 2 * W) - i, up * a[i]))
    win = a[best:bisect.bisect_right(a, a[best] + 2 * W)]; out['mode'] = win[len(win) // 2]
    near = lambda x: bisect.bisect_right(a, x + W) - bisect.bisect_left(a, x - W)
    out['sup2'] = next((x for x in s if near(x) >= 2), out['mode']); need = max(2, -(-n // 10))
    out['sup10'] = next((x for x in s if near(x) >= need), out['mode'])
    return out

def _scan(contig):
    import pysam
    want = _G['want'][contig]; out = collections.defaultdict(lambda: ([], []))
    for r in pysam.AlignmentFile(_G['bam']).fetch(contig):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or 'N' not in r.cigarstring: continue
        c = OL.read_chain(r)
        if c in want: o = out[c]; o[0].append(r.reference_start); o[1].append(r.reference_end)
    return contig, dict(out)

def bed_index(path, pad=0):
    d = collections.defaultdict(list)
    for l in open(path):
        f = l.split('\t')
        if len(f) > 2 and f[1].isdigit(): d[f[0]].append((int(f[1]) - pad, int(f[2]) + pad))
    return {c: (sorted(v), [a for a, _ in sorted(v)], max(b - a for a, b in v)) for c, v in d.items()}

def hit(ix, c, x):
    if c not in ix: return False
    v, st, mx = ix[c]; i = bisect.bisect_right(st, x) - 1
    while i >= 0 and st[i] >= x - mx:
        if v[i][1] > x: return True
        i -= 1
    return False

def main(gtf, bam, truth, out, *a):
    import pysam
    opt = dict(zip(a[::2], a[1::2])); thr = int(opt.get('--threads', 32))
    M = {}
    for c, s, ex, lines in oneread.models(gtf).values():
        if len(ex) > 1 and s in '+-' and 'one_read_' not in lines[0]: M[(c, oneread.chain(ex))] = s
    A = collections.defaultdict(lambda: ([], []))
    if truth != '-':
        for c, s, ex, _ in oneread.models(truth).values():
            if len(ex) > 1 and s in '+-': t = A[(c, s, oneread.chain(ex))]; t[0].append(ex[0][0] if s == '+' else ex[-1][1]); t[1].append(ex[-1][1] if s == '+' else ex[0][0])
    want = collections.defaultdict(set)
    for c, ch in M: want[c].add(ch)
    _G.update(bam=bam, want=want)
    with multiprocessing.get_context('fork').Pool(thr) as pool: R = dict(pool.map(_scan, [c for c in pysam.AlignmentFile(bam).references if c in want], chunksize=1))
    cage = bed_index(opt['--cage'], 25) if '--cage' in opt else None; pa = bed_index(opt['--polya'], 25) if '--polya' in opt else None
    K = ('ext', 'med', 'mode', 'sup2', 'sup10'); T = {e: {k: collections.Counter() for k in K} for e in '53'}; n = na = 0; moved = {e: {k: 0 for k in K} for e in '53'}
    with open(out, 'w') as o:
        o.write('contig\tstrand\treads\t' + '\t'.join(f'{e}_{k}' for e in '53' for k in K) + '\tannotated\n')
        for (c, ch), s in M.items():
            st, en = R.get(c, {}).get(ch, ((), ()))
            if len(st) < 2: continue
            n += 1; five = candidates(st, -1) if s == '+' else candidates(en, 1); three = candidates(en, 1) if s == '+' else candidates(st, -1)
            ann = A.get((c, s, ch)); na += ann is not None
            for e, cand, ix in (('5', five, cage), ('3', three, pa)):
                for k in K:
                    x = cand[k]; t = T[e][k]; moved[e][k] += abs(x - cand['ext']) > W
                    if ix is not None: t['signal'] += hit(ix, c, x); t['bg'] += hit(ix, c, x + 500)
                    if ann:
                        d = min(abs(x - y) for y in ann[e == '3']); t['a10'] += d <= 10; t['a50'] += d <= 50
            o.write(f'{c}\t{s}\t{len(st)}\t' + '\t'.join(str(five[k]) for k in K) + '\t' + '\t'.join(str(three[k]) for k in K) + f'\t{int(ann is not None)}\n')
    print(f'chains with >= 2 full-chain reads {n}; annotated {na}')
    print('end\tcandidate\tmoved > 10 bp from today\twithin 25 bp of a CAGE peak (5) / polyA site (3)\tthe same 500 bp away (background)\tannotated chains: within 10 bp\twithin 50 bp')
    for e in '53':
        for k in K: t = T[e][k]; print(e + "'", k, f'{moved[e][k] / n:.3f}', f"{t['signal'] / n:.3f}" if (cage if e == '5' else pa) is not None else 'NA', f"{t['bg'] / n:.3f}" if (cage if e == '5' else pa) is not None else 'NA', f"{t['a10'] / max(1, na):.3f}", f"{t['a50'] / max(1, na):.3f}", sep='\t')

def selftest():
    v = [100, 500, 502, 505, 505, 507, 900, 905]
    c = candidates(v, -1); assert c['ext'] == 100 and c['sup2'] == 500 and c['mode'] in (502, 505) and c['med'] == 505, c
    c = candidates(v, 1); assert c['ext'] == 905 and c['sup2'] == 905 and c['med'] == 505, c
    c = candidates([10] * 30 + [300], 1); assert c['sup10'] == 10 and c['sup2'] == 10 and c['ext'] == 300, c
    print('end_choice selftest ok')

if __name__ == '__main__': selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:])
