#!/usr/bin/env python3
"""PREREG row 1751: de novo local realignment of reads that skip an exon of <= 15 bp another read of the locus spells.

  skip_fix.py scan REF_GTF BAM FASTA OUT_PREFIX [THREADS] [R1746_TSV]   candidates, fixes, decoys, reach -> OUT.units.tsv, OUT.regions.bed
  skip_fix.py fix OUT.units.tsv BAM FASTA OUT_BAM                        rewrite the CIGARs -> OUT_BAM, OUT_BAM.fixes.tsv
  skip_fix.py score OUT.units.tsv NAME=GTF [...]                         emitted chains with the skip / the exon at fired units
  skip_fix.py --selftest
Candidate = a primary read's intron (d, a) where another primary read has the consecutive introns (d, x), (y, a) with
0 < y - x <= 15. A read is fixed when its bases around the junction fit the sequence with the exon better by MARGIN edits
(row 1746's test, window clamped to the read's own flanking exons). Only an exact (d, a) intron is a candidate: a read whose
junction the aligner also shifted is not reached."""
import sys, collections, multiprocessing, random, functools
from alt_events import load
from skip_realign import edit, W, MARGIN
MAXEX, DECOYS, DCAP, DMIN, FLANK = 15, 300, 20, 5, 100000

def introns(start, cig):
    pos, out = start, []
    for o, n in cig:
        if o == 3: out.append((pos, pos + n))
        if o in (0, 2, 3, 7, 8): pos += n
    return out, pos

def expand(cig):
    """one entry per base; an intron stays one (3, n) entry"""
    out = []
    for o, n in cig: out += [(3, n)] if o == 3 else [o] * n
    return out

def compress(ops):
    out = []
    for e in ops:
        if isinstance(e, tuple): out.append(e)
        elif out and out[-1][0] == e and e != 3: out[-1] = (e, out[-1][1] + 1)
        else: out.append((e, 1))
    return out

def cut(ops, start, lo, hi):
    """(i_lo, q_lo, i_hi, q_hi): entries [i_lo, i_hi) cover reference lo..hi; an insertion just before lo stays outside,
    one just before hi inside (qpos of skip_realign.py)"""
    r, q, a = start, 0, None
    for i, e in enumerate(ops):
        if a is None and r == lo and e in (0, 2, 7, 8): a = (i, q)
        if a is not None and r == hi and e != 1: return a + (i, q)
        if isinstance(e, tuple): r += e[1]
        else: r += e in (0, 2, 7, 8); q += e in (0, 1, 4, 7, 8)
    return a and a + (len(ops), q)

def align(q, t):
    """global edit alignment of read bases q to reference t: (distance, per-base ops 0 = M, 1 = I, 2 = D)"""
    n, m = len(q), len(t); D = [list(range(m + 1))] + [[i] + [0] * m for i in range(1, n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1): D[i][j] = min(D[i - 1][j] + 1, D[i][j - 1] + 1, D[i - 1][j - 1] + (q[i - 1] != t[j - 1]))
    ops, i, j = [], n, m
    while i or j:
        if i and j and D[i][j] == D[i - 1][j - 1] + (q[i - 1] != t[j - 1]): ops.append(0); i -= 1; j -= 1
        elif i and D[i][j] == D[i - 1][j] + 1: ops.append(1); i -= 1
        else: ops.append(2); j -= 1
    return D[n][m], ops[::-1]

def try_fix(cig, start, seq, cands, fetch):
    """the read's first candidate intron whose bases spell a short exon: (new cigar, (d, x, y, a)), else None"""
    ivs, end = introns(start, cig); ops = None
    for k, (d, a) in enumerate(ivs):
        if (d, a) not in cands: continue
        ops = ops or expand(cig); lo = max(ivs[k - 1][1] if k else start, d - W); hi = min(ivs[k + 1][0] if k + 1 < len(ivs) else end, a + W)
        c = cut(ops, start, lo, hi)
        if not c: continue
        seg = seq[c[1]:c[3]].upper(); up, dn = fetch(lo, d), fetch(a, hi)
        if seg == up + dn: continue                         # spells the skip exactly: cannot fit the exon better
        db, best = edit(seg, up + dn), None
        for x, y in cands[d, a]:
            da, al = align(seg, up + fetch(x, y) + dn)
            if da + MARGIN <= db and (best is None or da < best[0]): best = (da, al, x, y)
        if not best: continue
        _, al, x, y = best; mid, j = [], 0
        for e in al:
            mid.append(e)
            if e != 1:
                j += 1
                if j == len(up): mid.append((3, x - d))
                elif j == len(up) + y - x: mid.append((3, a - y))
        return compress(ops[:c[0]] + mid + ops[c[2]:]), (d, x, y, a)
    return None

def fix_all(cig, start, seq, cands, fetch):
    """every fixable intron of the read: (new cigar, [units]); the rewrite keeps the read's start, end and length"""
    done = []
    for _ in range(10):
        f = try_fix(cig, start, seq, cands, fetch)
        if not f: break
        assert introns(start, f[0])[1] == introns(start, cig)[1] and sum(n for o, n in f[0] if o in (0, 1, 4, 7, 8)) == len(seq)
        cig = f[0]; done.append(f[1])
    return cig, done

def primary_spliced(r): return not (r.is_secondary or r.is_supplementary or r.is_unmapped) and 'N' in r.cigarstring

def scan(job):
    bam, fa_p, contig = job; import pysam
    fa = pysam.FastaFile(fa_p); fetch = functools.lru_cache(100000)(lambda s, e: fa.fetch(contig, s, e).upper())
    n, pair = collections.Counter(), collections.Counter()
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if not primary_spliced(r): continue
        ivs = introns(r.reference_start, r.cigartuples)[0]; n.update(ivs)
        for (d, x), (y, a) in zip(ivs, ivs[1:]):
            if 0 < y - x <= MAXEX: pair[d, x, y, a] += 1
    cands = collections.defaultdict(list)
    for d, x, y, a in sorted(pair):
        if n[d, a]: cands[d, a].append((x, y))
    lens = sorted(y - x for v in cands.values() for x, y in v); rng = random.Random(f'1751 {contig}')
    plain = sorted(k for k, v in n.items() if v >= DMIN and k not in cands and k[1] - k[0] > 200); rng.shuffle(plain)
    # decoy: a plain intron given a made-up exon of a real candidate's length, cut from the middle of the intron
    decoy = {k: [((k[0] + k[1]) // 2, (k[0] + k[1]) // 2 + lens[i % len(lens)])] for i, k in enumerate(plain[:DECOYS])} if lens else {}
    fixed, dseen, dec = collections.Counter(), collections.Counter(), collections.Counter()
    if cands or decoy:
        for r in pysam.AlignmentFile(bam).fetch(contig):
            if not primary_spliced(r) or not r.query_sequence: continue
            cig, st, seq = r.cigartuples, r.reference_start, r.query_sequence
            if cands: fixed.update(fix_all(cig, st, seq, cands, fetch)[1])
            for iv in introns(st, cig)[0] if decoy else ():
                if iv in decoy and dseen[iv] < DCAP:
                    dseen[iv] += 1; L = decoy[iv][0][1] - decoy[iv][0][0]; dec[L, 'n'] += 1
                    dec[L, 'fired'] += bool(try_fix(cig, st, seq, {iv: decoy[iv]}, fetch))
    return [(contig, d, x, y, a, pair[d, x, y, a], n[d, a], fixed[d, x, y, a]) for (d, a), v in cands.items() for x, y in v], dec

def read_units(path):
    rows = [l.rstrip('\n').split('\t') for l in open(path)][1:]
    return [(r[0], *map(int, r[1:5]), r[5], *map(int, r[6:9])) for r in rows]   # contig d x y a class pair_reads skip_reads fixed

def selftest():
    rng = random.Random(1); ref = ''.join(rng.choice('ACGT') for _ in range(900)); fetch = lambda s, e: ref[s:e]
    assert compress(expand([(4, 2), (0, 5), (3, 90), (1, 1), (0, 3)])) == [(4, 2), (0, 5), (3, 90), (1, 1), (0, 3)]
    assert align('ACGT', 'ACGT') == (0, [0, 0, 0, 0]) and align('ACT', 'ACGT')[0] == 1 and align('', 'AC') == (2, [2, 2])
    cands = {(140, 500): [(300, 306)], (540, 800): [(650, 659)]}; seq = ref[100:140] + ref[300:306] + ref[500:540]
    cig, done = fix_all([(4, 0), (0, 40), (1, 6), (3, 360), (0, 40)][1:], 100, seq, cands, fetch)        # exon as an insertion
    assert done == [(140, 300, 306, 500)] and cig == [(0, 40), (3, 160), (0, 6), (3, 194), (0, 40)], cig
    cig, done = fix_all([(4, 3), (0, 40), (3, 360), (1, 6), (0, 40), (4, 2)], 100, 'TTT' + seq + 'GG', cands, fetch)   # clips kept
    assert cig == [(4, 3), (0, 40), (3, 160), (0, 6), (3, 194), (0, 40), (4, 2)], cig
    assert fix_all([(0, 40), (3, 360), (0, 40)], 100, ref[100:140] + ref[500:540], cands, fetch)[1] == []           # a real skip
    junk = ref[100:140] + ref[300:306].translate(str.maketrans('ACGT', 'TGCA')) + ref[500:540]   # six other bases at the junction
    assert fix_all([(0, 40), (1, 6), (3, 360), (0, 40)], 100, junk, cands, fetch)[1] == []
    two = ref[100:140] + ref[300:306] + ref[500:540] + ref[650:659] + ref[800:830]                                 # two in one read
    cig, done = fix_all([(0, 40), (1, 6), (3, 360), (0, 40), (3, 260), (1, 9), (0, 30)], 100, two, cands, fetch)
    assert len(done) == 2 and introns(100, cig)[0] == [(140, 300), (306, 500), (540, 650), (659, 800)], cig
    short = ref[130:140] + ref[300:306] + ref[500:540] + ref[800:830]                    # window clamped by the next intron
    cig, done = fix_all([(0, 10), (1, 6), (3, 360), (0, 40), (3, 260), (0, 30)], 130, short, cands, fetch)
    assert introns(130, cig)[0] == [(140, 300), (306, 500), (540, 800)], cig
    print('selftest ok')

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == '--selftest': selftest(); sys.exit()
    import pysam
    if mode == 'scan':
        ref_p, bam, fa, out = sys.argv[2:6]; thr = int(sys.argv[6]) if len(sys.argv) > 6 else 16; prior = sys.argv[7] if len(sys.argv) > 7 else None
        ann_i, ann_p = set(), set()
        for (c, s), d in load(ref_p).items():
            for ch in d:
                ann_i.update((c, *iv) for iv in ch); ann_p.update((c, d0, x, y, a) for (d0, x), (y, a) in zip(ch, ch[1:]))
        rows, dec = [], collections.Counter()
        with multiprocessing.Pool(thr) as pool:
            for r, d in pool.imap_unordered(scan, [(bam, fa, c) for c in pysam.AlignmentFile(bam).references]): rows += r; dec.update(d)
        klass = lambda c, d, x, y, a: {(1, 0): 'incl', (1, 1): 'both', (0, 1): 'skip', (0, 0): 'neither'}[(c, d, x, y, a) in ann_p, (c, d, a) in ann_i]
        rows = sorted((c, d, x, y, a, klass(c, d, x, y, a), p, n, f) for c, d, x, y, a, p, n, f in rows)
        with open(out + '.units.tsv', 'w') as o:
            o.write('contig\tdonor\texon_start\texon_end\tacceptor\tannotation\treads_with_exon\treads_with_skip\treads_fixed\n')
            for r in rows: o.write('\t'.join(map(str, r)) + '\n')
        with open(out + '.regions.bed', 'w') as o:
            cur = None
            for c, s, e in sorted((r[0], max(0, r[1] - FLANK), r[4] + FLANK) for r in rows if r[8]):
                if cur and cur[0] == c and s <= cur[2]: cur[2] = max(cur[2], e); continue
                if cur: o.write('%s\t%d\t%d\n' % tuple(cur))
                cur = [c, s, e]
            if cur: o.write('%s\t%d\t%d\n' % tuple(cur))
        print('annotation of the unit\tcandidate units\tfired units\tskip reads at fired units\treads fixed')
        for k in ('incl', 'both', 'skip', 'neither', 'all'):
            u = [r for r in rows if k in (r[5], 'all')]; f = [r for r in u if r[8]]
            print(f'{k}\t{len(u)}\t{len(f)}\t{sum(r[7] for r in f)}\t{sum(r[8] for r in f)}')
        n, f = sum(v for (L, k), v in dec.items() if k == 'n'), sum(v for (L, k), v in dec.items() if k == 'fired')
        print(f'decoy reads tested / fired (rate)\t{n}\t{f}\t{f / max(1, n):.4f}')
        print('decoy by exon length (tested:fired)\t' + ' '.join(f"{L}={dec[L, 'n']}:{dec[L, 'fired']}" for L in sorted({L for L, _ in dec})))
        if prior:
            want = {(r[0], int(r[1]), int(r[2])) for r in (l.split('\t') for l in list(open(prior))[1:]) if r[3] == '<=15' and sum(map(int, r[7:10])) and 2 * int(r[7]) >= sum(map(int, r[7:10]))}
            got = {(r[0], r[2], r[3]) for r in rows if r[8]}; cand = {(r[0], r[2], r[3]) for r in rows}
            print(f'row 1746 units (>= half the skip reads carry the exon)\t{len(want)}\tcandidates\t{len(want & cand)}\tfired\t{len(want & got)}\treach\t{len(want & got) / max(1, len(want)):.3f}')
    elif mode == 'fix':
        units, bam, fa_p, out = sys.argv[2:6]; fa = pysam.FastaFile(fa_p); cands = collections.defaultdict(lambda: collections.defaultdict(list))
        for c, d, x, y, a, *_ in read_units(units): cands[c][d, a].append((x, y))
        src = pysam.AlignmentFile(bam); n = 0
        fetch = functools.lru_cache(200000)(lambda c, s, e: fa.fetch(c, s, e).upper())
        with pysam.AlignmentFile(out, 'wb', template=src) as dst, open(out + '.fixes.tsv', 'w') as log:
            for r in src.fetch(until_eof=True):
                c = r.reference_name
                if c in cands and primary_spliced(r) and r.query_sequence:
                    cig, done = fix_all(r.cigartuples, r.reference_start, r.query_sequence, cands[c], lambda s, e, c=c: fetch(c, s, e))
                    if done:
                        log.write(f"{r.query_name}\t{c}\t{r.reference_start}\t{r.cigarstring}\t"); r.cigartuples = cig; n += 1
                        for t in ('NM', 'MD', 'cs'):       # stale after the rewrite (the assembler reads AS only)
                            if r.has_tag(t): r.set_tag(t, None)
                        log.write(f"{r.cigarstring}\t{';'.join('%d-%d' % (x, y) for _, x, y, _ in done)}\n")
                dst.write(r)
        print(f'reads rewritten\t{n}')
    elif mode == 'score':
        units = [u for u in read_units(sys.argv[2]) if u[8]]
        print('arm\tfired units, annotation has only the exon form: chains with the skip intron\tchains with the exon\tall fired units: with the skip intron\twith the exon')
        for name, p in (a.split('=', 1) for a in sys.argv[3:]):
            ivs, prs = collections.Counter(), collections.Counter()
            for (c, s), d in load(p).items():
                for ch in d:
                    ivs.update((c, *iv) for iv in ch); prs.update((c, d0, x, y, a) for (d0, x), (y, a) in zip(ch, ch[1:]))
            f = lambda us: f"{sum(ivs[c, d, a] for c, d, a in {(u[0], u[1], u[4]) for u in us})}\t{sum(prs[u[:5]] for u in us)}"
            print(f"{name}\t{f([u for u in units if u[5] == 'incl'])}\t{f(units)}")
