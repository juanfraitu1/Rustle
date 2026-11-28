#!/usr/bin/env python3
"""A second end for a chain: candidates and their judges (register 1993).

    alt_end_models.py GTF BAM GENOME.fa TRUTH.gtf|- OUT.tsv [--cage PEAKS.bed] [--polya SITES.bed] [--threads N]

See docs/PREREG_alt_end_models_2026-11-28.md for the rule. Prints one table; OUT.tsv = one line per secondary end."""
import sys, os, bisect, collections, multiprocessing, importlib.util
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread, end_choice as T
_s = importlib.util.spec_from_file_location('ec', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'end_choice.py')); EC = importlib.util.module_from_spec(_s); _s.loader.exec_module(EC)
W, INSIDE, MIN_READS, MIN_SHARE, PAD = 10, 50, 3, 0.20, 25
RC = str.maketrans('ACGTN', 'TGCAN')

def secondary(v, up, primary):
    """v = read end positions, up = +1 when a larger coordinate is more extended. -> (position, reads) of the best supported
    cluster at least INSIDE bp inside `primary`, or None"""
    c = collections.Counter(v); a = sorted(c); n = len(v); best = None
    for x in a:
        if up * (primary - x) < INSIDE: continue
        k = sum(c[y] for y in a[bisect.bisect_left(a, x - W):bisect.bisect_right(a, x + W)])
        if k >= MIN_READS and k >= MIN_SHARE * n and c[x] == max(c[y] for y in a[bisect.bisect_left(a, x - W):bisect.bisect_right(a, x + W)]):
            if best is None or (k, up * x) > (best[1], up * best[0]): best = (x, k)
    return best

def a_rich(fa, c, pos, plus):
    """the genome just past a 3' end, in transcript sense: >= 12 A in 20 bp or 6 in a row"""
    d = fa.fetch(c, pos, pos + 20).upper() if plus else fa.fetch(c, max(pos - 20, 0), pos).upper().translate(RC)[::-1]
    return d.count('A') >= 12 or 'AAAAAA' in d

def main(gtf, bam, fasta, truth, out, *a):
    import pysam
    opt = dict(zip(a[::2], a[1::2])); thr = int(opt.get('--threads', 32))
    M = {}
    for c, s, ex, lines in oneread.models(gtf).values():
        if len(ex) > 1 and s in '+-' and 'one_read_' not in lines[0]: M[(c, oneread.chain(ex))] = (s, ex)
    want = collections.defaultdict(set)
    for c, ch in M: want[c].add(ch)
    T._G.update(bam=bam, want=want)
    with multiprocessing.get_context('fork').Pool(thr) as pool: R = dict(pool.map(T._scan, [c for c in pysam.AlignmentFile(bam).references if c in want], chunksize=1))
    A = collections.defaultdict(lambda: ([], [])); any5 = collections.defaultdict(list); any3 = collections.defaultdict(list)
    if truth != '-':
        for c, s, ex, _ in oneread.models(truth).values():
            if s not in '+-': continue
            e5, e3 = (ex[0][0], ex[-1][1]) if s == '+' else (ex[-1][1], ex[0][0]); any5[(c, s)].append(e5); any3[(c, s)].append(e3)
            if len(ex) > 1: t = A[(c, s, oneread.chain(ex))]; t[0].append(e5); t[1].append(e3)
    for d in (any5, any3):
        for k in d: d[k].sort()
    near = lambda d, c, s, x: (lambda v, i: any(abs(v[j] - x) <= PAD for j in (i - 1, i) if 0 <= j < len(v)))(d.get((c, s), []), bisect.bisect_left(d.get((c, s), []), x))
    ix = {'5': EC.bed_index(opt['--cage'], PAD) if '--cage' in opt else None, '3': EC.bed_index(opt['--polya'], PAD) if '--polya' in opt else None}
    fa = pysam.FastaFile(fasta); S = collections.defaultdict(collections.Counter); n = 0
    with open(out, 'w') as o:
        o.write('contig\tstrand\tend\tprimary\tsecondary\treads\tchain_reads\ta_rich\tsignal\tannotated_site\n')
        for (c, ch), (s, ex) in M.items():
            st, en = R.get(c, {}).get(ch, ((), ()))
            if len(st) < 2: continue
            n += 1
            for end in '53':
                left = (end == '5') == (s == '+'); v, up = (st, -1) if left else (en, 1); p = T.supported(v, up); ann = any5 if end == '5' else any3
                k = S[end + ' primary']; k['n'] += 1; k['ann'] += near(ann, c, s, p)
                if ix[end] is not None: k['sig'] += EC.hit(ix[end], c, p); k['bg'] += EC.hit(ix[end], c, p + 500)
                sec = secondary(v, up, p)
                if not sec: continue
                ar = end == '3' and a_rich(fa, c, sec[0], s == '+'); cls = end + (' secondary, A-rich (removed)' if ar else ' secondary')
                k = S[cls]; k['n'] += 1; an = near(ann, c, s, sec[0]); k['ann'] += an; sg = EC.hit(ix[end], c, sec[0]) if ix[end] is not None else 'NA'
                if ix[end] is not None: k['sig'] += sg; k['bg'] += EC.hit(ix[end], c, sec[0] + 500)
                key = (c, s, ch)
                if key in A: k['chain_annotated'] += 1; k['own'] += any(abs(sec[0] - y) <= PAD for y in A[key][end == '3'])
                o.write(f'{c}\t{s}\t{end}\t{p}\t{sec[0]}\t{sec[1]}\t{len(v)}\t{int(ar)}\t{sg if sg == "NA" else int(sg)}\t{int(an)}\n')
    print(f'chains with >= 2 full-chain reads {n}')
    print('class\tn\tshare of the chains\tsignal (CAGE 5 / polyA 3)\tbackground\twithin 25 bp of an annotated transcript end\tannotated chains: within 25 bp of an end of the SAME chain')
    for k in sorted(S):
        t = S[k]; m = max(1, t['n']); sig = ix[k[0]] is not None
        print(k, t['n'], f"{t['n'] / n:.4f}", f"{t['sig'] / m:.3f}" if sig else 'NA', f"{t['bg'] / m:.3f}" if sig else 'NA', f"{t['ann'] / m:.3f}", f"{t['own']} of {t['chain_annotated']}" if 'secondary' in k else '', sep='\t')

def selftest():
    v = [1000] * 6 + [1002, 998] + [700] * 3 + [702] + [400]
    assert secondary(v, 1, 1000) == (700, 4) and secondary(v, 1, 720) is None and secondary([1000] * 20 + [700] * 3, 1, 1000) is None
    assert secondary([100] * 5 + [400] * 4, -1, 100) == (400, 4) and secondary([100] * 5 + [130] * 4, -1, 100) is None
    print('alt_end_models selftest ok')

if __name__ == '__main__': selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:])
