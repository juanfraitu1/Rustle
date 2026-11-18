#!/usr/bin/env python3
"""PREREG row 1714: annotated alternative-event pairs, expressed, and which tools emit both chains.

  alt_events.py REF_GTF BAM OUT_TSV NAME=GTF [NAME=GTF ...] [--altends GTF] [--threads N] [--contigs A,B]
  alt_events.py --selftest
Chains are coordinate-sorted 0-based half-open introns. Last stdout block is the class table."""
import sys, re, gzip, collections, itertools, multiprocessing
TID = re.compile(r'transcript_id "([^"]+)"')
op = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def load(path, want_attr=None):
    """(contig, strand) -> chain -> [(start, end, attr)] over its transcripts (multi-exon only)."""
    ex = collections.defaultdict(list); at = {}
    for l in op(path):
        f = l.split('\t')
        if len(f) < 9: continue
        m = TID.search(f[8])
        if not m: continue
        k = (f[0], f[6], m.group(1))
        if f[2] == 'exon': ex[k].append((int(f[3]) - 1, int(f[4])))
        elif f[2] == 'transcript' and want_attr:
            a = re.search(r'\b' + want_attr + r' "?(\d+)', f[8]); at[k] = int(a.group(1)) if a else None
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for k, e in ex.items():
        e.sort()
        if len(e) > 1: out[k[:2]][tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))].append((e[0][0], e[-1][1], at.get(k)))
    return out

def classify(chains, strand):
    """chains: chain -> (min start, max end). -> [(class, shift, chainA, chainB)], one row per unordered pair and class."""
    plus = strand == '+'; S = set(chains); out = set()
    slot = collections.defaultdict(list)
    for c in S:
        for i in range(len(c)): slot[(c[:i], c[i + 1:])].append(c[i])
    for (pre, suf), ins in slot.items():
        for x, y in itertools.combinations(sorted(ins), 2):
            if (x[0] == y[0]) == (x[1] == y[1]): continue
            a, b = pre + (x,) + suf, pre + (y,) + suf; n = len(a); i = len(pre)
            left = x[1] == y[1]                       # the lower-coordinate boundary moves
            donor = left == plus; t = i if plus else n - 1 - i   # intron index in transcript order
            sh = abs(x[0] - y[0]) if left else abs(x[1] - y[1])
            if donor and t > 0: k = 'SITE_int3'
            elif not donor and t < n - 1: k = 'SITE_int5'
            else:  # a terminal exon: its other end is the transcript end
                if left: ov = max(chains[a][0], chains[b][0]) < min(x[0], y[0])
                else: ov = min(chains[a][1], chains[b][1]) > max(x[1], y[1])
                k = ('SITE_first3' if ov else 'AFE') if donor else ('SITE_last5' if ov else 'ALE')
            out.add((k, sh, a, b))
    for c in S:
        for i in range(len(c) - 1):
            b = c[:i] + ((c[i][0], c[i + 1][1]),) + c[i + 2:]
            if b in S:
                out.add(('CASSETTE', 0, c, b))
                if c[i + 1][0] - c[i][1] <= 30: out.add(('MICRO', c[i + 1][0] - c[i][1], c, b))   # row 1744: skipped exon <= 30 bp
        for k in range(1, len(c)):
            if c[k:] in S: out.add(('NEST5' if plus else 'NEST3', 0, c, c[k:]))   # sub-chain reaches the high-coordinate end
            if c[:k] in S: out.add(('NEST3' if plus else 'NEST5', 0, c, c[:k]))
    for c in S:   # row 1744: retained intron = the same chain without one INTERNAL intron (a terminal one is NEST5 / NEST3)
        for i in range(1, len(c) - 1):
            if c[:i] + c[i + 1:] in S: out.add(('RI', c[i][1] - c[i][0], c, c[:i] + c[i + 1:]))
    mx = collections.defaultdict(list)   # row 1744: mutually exclusive = one internal exon replaced by a disjoint one, same outer sites
    for c in S:
        for i in range(len(c) - 1): mx[(c[:i], c[i][0], c[i + 1][1], c[i + 2:])].append((c[i][1], c[i + 1][0], c))
    for v in mx.values():
        for (_, a2, ca), (b1, _, cb) in itertools.combinations(sorted(v), 2):
            if a2 <= b1: out.add(('MXE', 0, ca, cb))
    return sorted(out)

def count(job):
    bam, contig = job; import pysam
    cnt = collections.Counter(); lo = collections.defaultdict(collections.Counter); hi = collections.defaultdict(collections.Counter)
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for o, ln in r.cigartuples:
            if o == 3: ch.append((pos, pos + ln))
            if o in (0, 2, 3, 7, 8): pos += ln
        ch = tuple(ch); cnt[ch] += 1; lo[ch][r.reference_start] += 1; hi[ch][pos] += 1
    return contig, cnt, lo, hi

def skip_bin(n): return '<=15' if n <= 15 else '16-30' if n <= 30 else '31-60' if n <= 60 else '>60'

def selftest():
    A = ((100, 200), (300, 400), (500, 600)); ends = lambda *cs: {c: (0, 1000) for c in cs}
    g = lambda cs, s: {(k, a, b) for k, _, a, b in classify(ends(*cs), s)}
    B = ((100, 200), (310, 400), (500, 600))    # middle intron: low boundary moves
    assert g((A, B), '+') == {('SITE_int3', A, B)} and g((A, B), '-') == {('SITE_int5', A, B)}
    C = ((100, 200), (300, 400), (500, 620))    # last intron in coordinates, high boundary moves
    assert g((A, C), '+') == {('SITE_last5', A, C)} and g((A, C), '-') == {('SITE_first3', A, C)}
    assert {k for k, _, _, _ in classify({A: (0, 1000), C: (0, 610)}, '-')} == {'AFE'}   # first exons (610.., 620..) do not overlap
    D = ((100, 400), (500, 600))
    assert g((A, D), '+') == {('CASSETTE', A, D)}
    assert g((A, A[1:]), '+') == {('NEST5', A, A[1:])} and g((A, A[1:]), '-') == {('NEST3', A, A[1:])}
    assert classify(ends(A, B), '+')[0][1] == 10
    E4 = ((100, 200), (300, 400), (500, 600), (700, 800))
    assert g((E4, E4[:1] + E4[2:]), '+') == {('RI', E4, E4[:1] + E4[2:])}   # internal intron (300, 400) retained
    M1 = ((100, 200), (220, 400), (500, 600)); M2 = ((100, 300), (320, 400), (500, 600)); SK = ((100, 400), (500, 600))
    assert g((M1, M2), '+') == {('MXE', M1, M2)}   # exons (200, 220) and (300, 320), same outer sites 100 and 400
    assert g((M1, SK), '+') == {('CASSETTE', M1, SK), ('MICRO', M1, SK)}   # 20 bp exon skipped
    assert skip_bin(15) == '<=15' and skip_bin(30) == '16-30' and skip_bin(60) == '31-60' and skip_bin(61) == '>60'
    print('selftest ok')

if __name__ == '__main__':
    if sys.argv[1] == '--selftest': selftest(); sys.exit()
    a = sys.argv[1:]; opt = lambda k, d=None: a[a.index(k) + 1] if k in a else d
    ref_p, bam, out_p = a[:3]; arms = [x.split('=', 1) for x in a[3:] if '=' in x and not x.startswith('--')]
    thr = int(opt('--threads', 16)); only = set(opt('--contigs', '').split(',')) - {''}
    import pysam
    ref = load(ref_p); have = set(pysam.AlignmentFile(bam).references)
    contigs = sorted({c for c, _ in ref if c in have and (not only or c in only)})
    emitted = {}
    for n, p in arms:
        e = collections.defaultdict(set)
        for (c, _), d in load(p).items(): e[c] |= set(d)
        emitted[n] = e
    alt = load(opt('--altends'), 'alt_tes') if opt('--altends') else None
    with multiprocessing.Pool(thr) as pool: reads = {c: x for c, *x in pool.imap_unordered(count, [(bam, c) for c in contigs])}
    near = lambda ctr, p, w=25: sum(n for q, n in ctr.items() if abs(q - p) <= w)
    T = collections.defaultdict(collections.Counter); same = collections.Counter(); SKIP = collections.defaultdict(collections.Counter)
    with open(out_p, 'w') as o:
        o.write('contig\tstrand\tclass\tshift\treads_a\treads_b\t' + '\t'.join(n for n, _ in arms) + '\tchain_a\tchain_b\n')
        for (c, s), d in sorted(ref.items()):
            if c not in reads: continue
            cnt, lo, hi = reads[c]
            for k, sh, x, y in classify({ch: (min(t[0] for t in v), max(t[1] for t in v)) for ch, v in d.items()}, s):
                rx, ry = cnt.get(x, 0), cnt.get(y, 0)
                if rx < ry: x, y, rx, ry = y, x, ry, rx   # x = major
                if k.startswith('SITE'): k += '<=10' if sh <= 10 else '>10'
                T[k]['annotated'] += 1
                if ry < 2: continue
                st = [('all', True), ('minor2-4', ry < 5), ('minor>=5', ry >= 5), ('share<0.1', ry < 0.1 * (rx + ry)), ('share>=0.1', ry >= 0.1 * (rx + ry))]
                em = []
                for n, _ in arms:
                    ex_, ey = x in emitted[n][c], y in emitted[n][c]; v = 'both' if ex_ and ey else 'major' if ex_ else 'minor' if ey else 'neither'
                    em.append(v)
                    for sn, ok in st:
                        if ok: T[k][(sn, n, v)] += 1
                for sn, ok in st:
                    if ok: T[k][(sn, 'n')] += 1
                o.write(f'{c}\t{s}\t{k}\t{sh}\t{rx}\t{ry}\t' + '\t'.join(em) + f'\t{x}\t{y}\n')
            for ch in d:   # row 1744: reads that skip an internal exon no annotated chain skips (unit = chain x exon, >= 5 reads)
                rc = cnt.get(ch, 0)
                for i in range(len(ch) - 1):
                    b = ch[:i] + ((ch[i][0], ch[i + 1][1]),) + ch[i + 2:]
                    if b in d: continue
                    rb = cnt.get(b, 0)
                    if rc + rb < 5: continue
                    K = SKIP[skip_bin(ch[i + 1][0] - ch[i][1])]; K['units'] += 1; K['full'] += rc; K['skip'] += rb; K['mostly'] += rb >= rc
                    for n, _ in arms:
                        K[n, 'full'] += ch in emitted[n][c]; K[n, 'skip_only'] += b in emitted[n][c] and ch not in emitted[n][c]
            for ch, v in d.items():   # SAME3 / SAME5: one chain, annotated ends >= 100 bp apart
                for side in ('3', '5'):
                    low = (side == '5') == (s == '+'); ann = sorted({t[0] if low else t[1] for t in v}); ctr = (lo if low else hi).get(ch, {})
                    if ann[-1] - ann[0] < 100: continue
                    same[(side, 'annotated')] += 1
                    p, q = ann[0], ann[-1]
                    if near(ctr, p) < 3 or near(ctr, q) < 3: continue
                    same[(side, 'expressed')] += 1
                    if side == '3' and alt is not None:
                        m = alt.get((c, s), {}).get(ch) or [t for (c2, _), dd in alt.items() if c2 == c for t in dd.get(ch, [])]
                        if not m: same[('3', 'chain not emitted')] += 1; continue
                        e3 = m[0][1] if s == '+' else m[0][0]; a3 = m[0][2]
                        a3 = None if a3 is None else (a3 if s == '+' else a3 - 1)   # alt_tes is 1-based
                        hit = lambda z, w: z is not None and abs(z - w) <= 25
                        same[('3', 'emitted end at one annotated end')] += hit(e3, p) or hit(e3, q)
                        same[('3', 'emitted at one, alt_tes at the other')] += (hit(e3, p) and hit(a3, q)) or (hit(e3, q) and hit(a3, p))
    names = [n for n, _ in arms]
    for sn in ('all', 'minor2-4', 'minor>=5', 'share<0.1', 'share>=0.1'):
        print(f'## stratum {sn}: class, annotated pairs, expressed pairs, then per arm both/major/minor/neither (both share)')
        for k in sorted(T):
            n = T[k][(sn, 'n')]
            print(f"{k}\t{T[k]['annotated']}\t{n}\t" + '\t'.join(f"{m} {'/'.join(str(T[k][(sn, m, v)]) for v in ('both', 'major', 'minor', 'neither'))} ({T[k][(sn, m, 'both')] / max(1, n):.3f})" for m in names))
    print('## reads skipping an internal exon that no annotated chain skips: exon length, chain x exon units, full-chain reads, skip reads (share), units with skip >= full, then per arm units with the full chain emitted / only the skip chain emitted')
    for k in ('<=15', '16-30', '31-60', '>60'):
        K = SKIP[k]; print(f"SKIP{k}\t{K['units']}\t{K['full']}\t{K['skip']} ({K['skip'] / max(1, K['full'] + K['skip']):.4f})\t{K['mostly']}\t" + '\t'.join(f"{n} {K[n, 'full']}/{K[n, 'skip_only']}" for n in names))
    print('## same chain, two annotated ends >= 100 bp apart'); [print(f'SAME{k[0]}\t{k[1]}\t{v}') for k, v in sorted(same.items())]
