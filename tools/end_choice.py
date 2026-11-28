#!/usr/bin/env python3
"""Supported ends (register 1992; driver: RUSTLE_SUPPORTED_ENDS in `all`).

    end_choice.py IN.gtf BAM OUT.gtf [0|3|5|both] [THREADS] [GENOME.fa]   |   end_choice.py --selftest

The assembler ends a transcript at the most extended end among the reads with exactly its intron chain: one stray read sets
it. Here each chosen end of a multi-exon transcript moves to the most extended read end that at least max(2, 10% of the
chain's reads) reach within W bp (reads = primary alignments with exactly the chain; none such: the densest 2 W window).
An end only moves INWARD and stays inside its exon; transcripts with fewer than two such reads, single-exon transcripts
and one-read transcripts (one_read_chain / one_read_locus) are untouched. Moved transcript lines get ends_moved "5", "3" or "53".
Why (row 1992, held-out libraries): 3' end within 10 bp of the annotated end 0.54 -> 0.57 / 0.52 -> 0.59 (mouse adrenal,
cortex), 0.49 -> 0.68 (gorilla KB3781), at a polyA site 0.75 -> 0.80 (WTC11 r2); 5' end near a CAGE peak 0.76 -> 0.84 and
within 50 bp of the annotated start 0.66 -> 0.81 (WTC11 r2), 0.65 -> 0.69 / 0.56 -> 0.66 (mouse). Costs: simulated reads
0.005 (3') and 0.02-0.03 (5') at 10 bp; gorilla RefSeq starts 0.61 -> 0.47 (predicted models at the most upstream read).
With GENOME.fa (row 1993; driver: RUSTLE_ALT_ENDS=1) a chain also gets a SECOND model, ending at its main 3' end, when a
cluster of read ends (within W bp of its commonest position) lies >= 50 bp inside the first model's 3' end, holds >= 3 reads
and >= half of the chain's reads, has a polyA signal (AATAAA / ATTAAA) 5-40 bp upstream and no A-rich genome after it
(>= 12 A in 20 bp or 6 in a row = internal priming). Its id ends in .a3, its line carries alt_end "3". Held out: such ends
are at a polyA site 0.70 (WTC11 r2) and 0.66 (HG002) against 0.80 / 0.46 for first ends; 3-5% of the chains; under 0.2% on
simulated reads. The looser version (>= 20% of the reads, no signal) was 0.51-0.61 and is not offered."""
import bisect, collections, multiprocessing, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import oneread, one_read_loci as OL

W = 10
_G = {}

def supported(v, up):
    """v = read end positions; up = +1 when a larger coordinate is the more extended end, -1 otherwise"""
    a = sorted(v); n = len(a); need = max(2, -(-n // 10))
    for x in (reversed(a) if up > 0 else a):
        if bisect.bisect_right(a, x + W) - bisect.bisect_left(a, x - W) >= need: return x
    best = max(range(n), key=lambda i: (bisect.bisect_right(a, a[i] + 2 * W) - i, up * a[i]))
    win = a[best:bisect.bisect_right(a, a[best] + 2 * W)]; return win[len(win) // 2]

INSIDE, ALT_READS, ALT_SHARE = 50, 3, 0.5
RC = str.maketrans('ACGTN', 'TGCAN')

def secondary(v, up, primary):
    """(position, reads) of the best supported cluster of read ends >= INSIDE bp inside `primary`, or None"""
    c = collections.Counter(v); a = sorted(c); n = len(v); best = None
    for x in a:
        if up * (primary - x) < INSIDE: continue
        near = a[bisect.bisect_left(a, x - W):bisect.bisect_right(a, x + W)]; k = sum(c[y] for y in near)
        if k >= ALT_READS and k >= ALT_SHARE * n and c[x] == max(c[y] for y in near) and (best is None or (k, up * x) > (best[1], up * best[0])): best = (x, k)
    return best

def polya_context(fetch, c, pos, plus):
    """(polyA signal 5-40 bp upstream, A-rich genome after) of a 3' end, in transcript sense; fetch(contig, start, end)"""
    if plus: up, down = fetch(c, max(pos - 40, 0), max(pos - 5, 0)).upper(), fetch(c, pos, pos + 20).upper()
    else: up, down = fetch(c, pos + 5, pos + 40).upper().translate(RC)[::-1], fetch(c, max(pos - 20, 0), pos).upper().translate(RC)[::-1]
    return 'AATAAA' in up or 'ATTAAA' in up, down.count('A') >= 12 or 'AAAAAA' in down

def _scan(contig):
    import pysam
    want = _G['want'][contig]; out = collections.defaultdict(lambda: ([], []))
    for r in pysam.AlignmentFile(_G['bam']).fetch(contig):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or 'N' not in r.cigarstring: continue
        c = OL.read_chain(r)
        if c in want: o = out[c]; o[0].append(r.reference_start); o[1].append(r.reference_end)
    return contig, dict(out)

def new_ends(ex, strand, starts, ends, which):
    """-> (start, end) 0-based half-open of the transcript after the rule; inward only, inside the end exons"""
    a, b = ex[0][0], ex[-1][1]
    if len(starts) < 2: return a, b
    left, right = ('5', '3') if strand == '+' else ('3', '5')
    if left in which: a = min(max(a, supported(starts, -1)), ex[0][1] - 1)
    if right in which: b = max(min(b, supported(ends, 1)), ex[-1][0] + 1)
    return a, b

def main(gtf, bam, out, which='both', threads=4, fasta=''):
    import pysam
    which = {'both': '53', '3': '3', '5': '5', '0': ''}[which]; fa = pysam.FastaFile(fasta) if fasta else None; alt = {}
    M = oneread.models(gtf); want = collections.defaultdict(set)
    for c, s, ex, lines in M.values():
        if len(ex) > 1 and s in '+-' and 'one_read_' not in lines[0]: want[c].add(oneread.chain(ex))
    _G.update(bam=bam, want=want)
    contigs = [c for c in pysam.AlignmentFile(bam).references if c in want]
    with multiprocessing.get_context('fork').Pool(max(1, min(int(threads), len(contigs) or 1))) as pool: R = dict(pool.map(_scan, contigs, chunksize=1))
    move = {}
    for t, (c, s, ex, lines) in M.items():
        if len(ex) < 2 or s not in '+-' or 'one_read_' in lines[0]: continue
        st, en = R.get(c, {}).get(oneread.chain(ex), ((), ())); a, b = new_ends(ex, s, st, en, which)
        if fa is not None and len(st) >= 2:
            v, up, prim = (en, 1, b) if s == '+' else (st, -1, a); sec = secondary(v, up, prim)
            if sec:
                sig, rich = polya_context(fa.fetch, c, sec[0], s == '+')
                if sig and not rich: alt[t] = (sec[0], sec[1], len(v))
        if (a, b) != (ex[0][0], ex[-1][1]):
            left, right = a != ex[0][0], b != ex[-1][1]       # which genomic side moved; the strand says which end that is
            move[t] = (ex[0], ex[-1], a, b, ''.join(e for e, hit in (('5', left if s == '+' else right), ('3', right if s == '+' else left)) if hit))
    n5 = sum('5' in m[4] for m in move.values()); n3 = sum('3' in m[4] for m in move.values())
    extra = collections.defaultdict(list)       # the second model of a chain: its lines, written after the first model's
    with open(out + '.tmp', 'w') as o:
        def flush(t):
            for x in extra.pop(t, ()): o.write(x)
        last_t = None
        for l in open(gtf):
            f = l.split('\t'); m = re.search(r'transcript_id "([^"]+)"', f[8]) if len(f) > 8 else None; mv = move.get(m.group(1)) if m else None
            t = m.group(1) if m else None
            if t != last_t and last_t is not None: flush(last_t)
            last_t = t
            if t in alt and f[2] in ('transcript', 'exon'):
                c, s, ex, _ = M[t]; pos, k, tot = alt[t]; g = f[:]; g[8] = g[8].replace(f'transcript_id "{t}"', f'transcript_id "{t}.a3"')
                lastx = ex[-1] if s == '+' else ex[0]; e = (int(g[3]) - 1, int(g[4]))
                if g[2] == 'transcript' or e == lastx:
                    if s == '+': g[4] = str(pos)
                    else: g[3] = str(pos + 1)
                if g[2] == 'transcript': g[8] = g[8].rstrip('\n').rstrip(';') + f'; alt_end "3"; alt_end_reads "{k}/{tot}";\n'
                extra[t].append('\t'.join(g))
            if mv:
                first, last, a, b, tag = mv
                if f[2] == 'transcript': f[3], f[4] = str(a + 1), str(b); f[8] = f[8].rstrip('\n').rstrip(';') + f'; ends_moved "{tag}";\n'
                elif f[2] == 'exon':
                    e = (int(f[3]) - 1, int(f[4]))
                    if e == first: f[3] = str(a + 1)
                    if e == last: f[4] = str(b)
                l = '\t'.join(f)
            o.write(l)
        if last_t is not None: flush(last_t)
    os.replace(out + '.tmp', out)
    print(f"supported ends: {len(move)} of {sum(len(v) for v in want.values())} spliced chains moved (5' end {n5}, 3' end {n3}; W = {W} bp)" + (f"; second models for a main 3' end {len(alt)}" if fa is not None else ''))

def selftest():
    assert supported([100, 500, 502, 505, 505, 507, 900, 905], -1) == 500 and supported([100, 500, 502, 505, 505, 507, 900, 905], 1) == 905
    assert supported([10] * 30 + [300], 1) == 10 and supported([5, 400, 900], -1) == 5          # no supported end: densest window, most extended on ties
    ex = [(100, 200), (300, 400), (500, 900)]
    assert new_ends(ex, '+', [100, 150, 152], [900, 700, 705], '53') == (150, 705) and new_ends(ex, '+', [100, 150, 152], [900, 700, 705], '3') == (100, 705)
    assert new_ends(ex, '-', [100, 150, 152], [900, 700, 705], '5') == (100, 705) and new_ends(ex, '+', [100], [900], '53') == (100, 900)
    assert new_ends(ex, '+', [50, 60, 62], [950, 940, 941], '53') == (100, 900)                 # never outward
    v = [1000] * 4 + [700] * 5 + [702]
    assert secondary(v, 1, 1000) == (700, 6) and secondary([1000] * 20 + [700] * 5, 1, 1000) is None and secondary([100] * 4 + [400] * 5, -1, 100) == (400, 5)
    seq = 'C' * 60 + 'AATAAA' + 'C' * 14 + 'G' * 40; f = lambda c, a, b: seq[a:b]
    assert polya_context(f, 'c', 80, True) == (True, False) and polya_context(lambda c, a, b: ('C' * 80 + 'A' * 20)[a:b], 'c', 80, True) == (False, True)
    print('end_choice selftest ok')

if __name__ == '__main__': selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:])
