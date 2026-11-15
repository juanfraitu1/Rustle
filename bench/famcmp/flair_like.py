#!/usr/bin/env python3
"""PREREG r1640 (docs/PREREG_flair_profile_2026-11-11.md): chain-level emulation of FLAIR collapse, without realignment.

  flair_like.py BAM REFSEQ DEFAULT_GTF FLAIR_GTF OUT_TSV [THREADS]      (selftest: flair_like.py selftest)
OUT_TSV: contig, introns, reads with the chain, assigned reads (E1), dropped as subset, assigned + unspliced reads (E2), in refseq, in flair."""
import sys, collections, multiprocessing, statistics
sys.path.insert(0, __import__('os').path.dirname(__file__))
from tool_match import read
TOL, WIN, MIN = 10, 100, 3
FP = False  # workers return FLAIR's first pass (E4) instead of the chain-level arms

def collapse(reads, mono=()):
    """reads: [(start, end, chain)] of one contig, mono: [(start, end)] unspliced reads
    -> {chain: [reads with the chain, assigned reads (E1), dropped as subset, assigned + unspliced reads (E2)]}."""
    g = collections.defaultdict(lambda: ([], []))
    for a, b, ch in reads: g[ch][0].append(a); g[ch][1].append(b)
    info = {ch: (len(s), int(statistics.median(s)), int(statistics.median(e))) for ch, (s, e) in g.items()}
    by = collections.defaultdict(set)
    for ch in info:
        for i in ch: by[i].add(ch)
    def frame(c, s):  # the exon ends of c around the run s inside c, or None
        n = len(s)
        for k in range(len(c) - n + 1):
            if c[k:k + n] == s: return (c[k - 1][1] if k else info[c][1], c[k + n][0] if k + n < len(c) else info[c][2])
    drop = set()
    for s, (n, a, b) in info.items():
        sup = []
        for c in by[s[0]]:
            if len(c) <= len(s): continue
            f = frame(c, s)
            if f and (a > f[0] - WIN or b < f[1] + WIN): sup.append(info[c][0])
        if sup and (len(s) < 3 or not (n > 3 and n > 1.2 * max(sup))): drop.add(s)
    tlen = lambda c: info[c][2] - info[c][1] - sum(y - x for x, y in c)
    out = {ch: [v[0], 0, int(ch in drop)] for ch, v in info.items()}
    for a, b, s in reads:
        if s not in drop: out[s][1] += 1; continue
        best = None
        for c in by[s[0]]:
            if c in drop or len(c) <= len(s): continue
            f = frame(c, s)
            if f and a >= f[0] - TOL and b <= f[1] + TOL and (best is None or (tlen(c), c) < best): best = (tlen(c), c)
        if best: out[best[1]][1] += 1
    for c, v in out.items(): v.append(v[1])
    bins = collections.defaultdict(list)  # 1 kb bins of the kept chains' exons: (start, end, transcript length, chain)
    for c, (n, a, b) in info.items():
        if c in drop: continue
        t = tlen(c)
        for x, y in zip((a,) + tuple(i[1] for i in c), tuple(i[0] for i in c) + (b,)):
            for k in range(x // 1000, y // 1000 + 1): bins[k].append((x, y, t, c))
    for a, b in mono:
        best = min(((t, c) for x, y, t, c in bins.get(a // 1000, ()) if x <= a + TOL and y >= b - TOL), default=None)
        if best: out[best[1]][3] += 1
    for c, v in out.items(): v += [info[c][1], info[c][2]]  # median ends, for the first-pass sequences (E3)
    return out

def best_sites(sites, total, tss):
    """FLAIR's find_tsss: up to 2 end sites of one chain. sites {position: reads}; -> [(position, reads exactly there)]."""
    sites = dict(sites); found = []; remaining = sum(sites.values()); avg = remaining / len(sites)
    while sites and remaining / total > 0.25 and len(found) < 2:
        best = None
        for s in sorted(sites, reverse=not tss):  # ties to the longer end
            w = sum(n if s == t else (WIN - abs(s - t)) / WIN * n for t, n in sites.items() if abs(s - t) < WIN)
            if best is None or w > best[0]: best = (w, s, sites[s])
        left = {t: n for t, n in sites.items() if abs(t - best[1]) > WIN}; used = remaining - sum(left.values())
        sites, remaining = left, remaining - used
        if found and (used / total < 0.25 or best[2] < 4 or used < avg): break
        found.append((best[1], best[2]))
    return found

def first_pass(reads, mono=()):
    """FLAIR's first pass (r1640 E4/E5): [(chain, start, end)] kept after collapse_isoforms_precise's end choice and the default
    subset filter; with mono (unspliced reads) also its single-exon isoforms, chain ()."""
    g = collections.defaultdict(list)
    for a, b, ch in reads: g[ch].append((a, b))
    iso = collections.defaultdict(list)  # chain -> [(start, end, support)]
    for ch, se in g.items():
        tss = collections.Counter(a for a, b in se)
        for a, na in best_sites(tss, len(se), True):
            tes = collections.Counter(b for x, b in se if abs(x - a) <= WIN)
            for b, nb in best_sites(tes, len(se), False): iso[ch].append((a, b, max(na, nb)))
    by = collections.defaultdict(set)
    for ch in iso:
        for i in ch: by[i].add(ch)
    keep = []
    for s, ends in iso.items():
        big = [(c, k) for c in by[s[0]] if len(c) > len(s) for k in range(len(c) - len(s) + 1) if c[k:k + len(s)] == s]
        for a, b, n in ends:
            sup = []
            for c, k in big:
                for ca, cb, cn in iso[c]:
                    lo = c[k - 1][1] if k else ca; hi = c[k + len(s)][0] if k + len(s) < len(c) else cb
                    sup += [cn] * ((a > lo - WIN) + (b < hi + WIN))
            if not sup or (len(s) >= 3 and n > 3 and n > 1.2 * max(sup)): keep.append((s, a, b))
    if not mono: return keep
    groups = []  # [mean start, mean end, reads, [(start, end)]]: unspliced reads joined while they start inside the group
    for a, b in sorted(mono):
        g = groups[-1] if groups else None
        if g and a <= g[1] and b <= g[1] + WIN: g[2] += 1; g[0] += (a - g[0]) / g[2]; g[1] += (b - g[1]) / g[2]; g[3].append((a, b))
        else: groups.append([float(a), float(b), 1, [(a, b)]])
    se = []
    for _, _, n, rs in groups:
        t = best_sites(collections.Counter(a for a, b in rs), n, True)[:1]
        if not t: continue
        e = best_sites(collections.Counter(b for a, b in rs if abs(a - t[0][0]) <= WIN), n, False)[:1]
        if e and e[0][0] > t[0][0]: se.append((t[0][0], e[0][0], max(t[0][1], e[0][1])))
    bins = collections.defaultdict(set); sbin = collections.defaultdict(list)
    for ch, ends in iso.items():
        for a, b, _ in ends:
            for x, y in zip((a,) + tuple(i[1] for i in ch), tuple(i[0] for i in ch) + (b,)):
                for k in range((x - WIN) // 1000, (y + WIN) // 1000 + 1): bins[k].add((x, y))
    for a, b, n in se:
        for k in range((a - WIN) // 1000, (b + WIN) // 1000 + 1): sbin[k].append((a, b, n))
    for a, b, n in se:
        if any(x - WIN <= a and b <= y + WIN for x, y in bins.get(a // 1000, ())): continue
        big = [m for x, y, m in sbin.get(a // 1000, ()) if y - x > b - a and x - WIN <= a and b <= y + WIN]
        if not big or (n > 3 and all(m * 1.2 < n for m in big)): keep.append(((), a, b))
    return keep

def assign(sam, per_read=None):
    """FLAIR's read assignment (count_sam_transcripts.getbesttranscript, no --stringent): per read, drop alignments with an
    indel > 25 bp, then best AS, least clipping, shortest transcript, name. -> (reads per transcript, reads, assigned)."""
    import re
    tl = {}; cnt = collections.Counter(); cur = None; best = None; reads = assigned = 0
    for l in sam:
        if l[0] == '@':
            if l.startswith('@SQ'): f = dict(x.split(':', 1) for x in l.rstrip().split('\t')[1:]); tl[f['SN']] = int(f['LN'])
            continue
        f = l.split('\t', 12); flag = int(f[1])
        if f[0] != cur:
            reads += 1
            if best:
                cnt[best[3]] += 1; assigned += 1
                if per_read: per_read.write(f'{cur}\t{best[3]}\n')
            cur, best = f[0], None
        if flag & 0x804 or f[2] == '*': continue
        ops = re.findall(r'(\d+)([MIDNSHP=X])', f[5])
        if any(o in 'IDN' and int(n) > 25 for n, o in ops): continue
        a = re.search(r'\tAS:i:(-?\d+)', l); k = (-int(a.group(1)), sum(int(n) for n, o in ops if o in 'SH'), tl[f[2]], f[2])
        if best is None or k < best: best = k
    if best:
        cnt[best[3]] += 1; assigned += 1
        if per_read: per_read.write(f'{cur}\t{best[3]}\n')
    return cnt, reads, assigned

def gtf(ids, cnt, w):
    """ids {id: (contig, strand, start, end, chain)}, cnt {id: assigned reads} -> GTF of the isoforms with >= MIN reads."""
    for i, (c, st, a, b, ch) in ids.items():
        n = cnt.get(i, 0)
        if n < MIN: continue
        at = f'gene_id "FL.{i}"; transcript_id "FL.{i}"; reads "{n}";'
        w.write(f'{c}\tflairlike\ttranscript\t{a + 1}\t{b}\t.\t{st}\t.\t{at}\n')
        for x, y in zip((a,) + tuple(k[1] for k in ch), tuple(k[0] for k in ch) + (b,)): w.write(f'{c}\tflairlike\texon\t{x + 1}\t{y}\t.\t{st}\t.\t{at}\n')

def work(job):
    import pysam
    bam, contig = job; reads = []
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for op, ln in r.cigartuples:
            if op == 3: ch.append((pos, pos + ln))
            if op in (0, 2, 3, 7, 8): pos += ln
        reads.append((r.reference_start, pos, tuple(ch)))
    mono = [(r.reference_start, r.reference_end) for r in pysam.AlignmentFile(bam).fetch(contig) if not (r.is_secondary or r.is_supplementary or r.is_unmapped) and 'N' not in r.cigarstring]
    if FP: return contig, first_pass(reads, mono if FP == 2 else ())
    return contig, collapse(reads, mono)

if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        A = ((100, 200), (300, 400), (500, 600), (700, 800)); S = ((300, 400), (500, 600)); L = ((300, 400), (500, 600), (700, 800))
        r = collapse([(0, 900, A)] + [(250, 650, S)] * 2 + [(250, 900, L)] * 5)
        assert r[S][:3] == [2, 0, 1] and r[A][1] == 1, r   # S (3 exons) is a subset: dropped
        assert r[L][:3] == [5, 7, 0], r                    # L: 4 exons, 5 reads > 3 and > 1.2 x 1: kept, and the shortest container of S's reads
        r = collapse([(0, 900, A)] * 5 + [(250, 900, L)] * 5); assert r[L][2] == 1 and r[A][1] == 10, r  # 5 is not > 1.2 x 5
        r = collapse([(0, 900, A)] + [(250, 900, L)] * 5, [(520, 590), (0, 90), (650, 690)]); assert r[A][3] == 2 and r[L][3] == 6 and r[L][1] == 5, r  # unspliced reads: shortest chain with a containing exon
        sam = ['@SQ\tSN:1\tLN:900\n', '@SQ\tSN:2\tLN:500\n'] + [f'{q}\t{fl}\t{t}\t1\t0\t{cg}\t*\t0\t0\t*\t*\tAS:i:{sc}\n' for q, fl, t, cg, sc in
               (('r1', 0, '1', '100M', 90), ('r1', 256, '2', '100M', 90), ('r2', 0, '1', '50M30D50M', 99), ('r2', 256, '2', '100M', 80), ('r3', 4, '*', '*', 0), ('r4', 0, '1', '100M', 70), ('r4', 256, '2', '90M10S', 70))]
        assert best_sites({100: 5, 120: 1, 900: 4}, 10, True) == [(100, 5), (900, 4)] and best_sites({100: 5, 900: 3}, 8, True) == [(100, 5)]  # a second site needs 4 reads on it
        fp = first_pass([(0, 900, A)] * 5 + [(250, 900, L)] * 5 + [(250, 650, S)] * 9); assert sorted(c for c, a, b in fp) == [A], fp  # L: 5 is not > 1.2 x 5; S: < 4 exons
        assert sorted(c for c, a, b in first_pass([(0, 900, A)] + [(250, 900, L)] * 5)) == [A, L]
        fp = first_pass([(0, 900, A)], [(520, 590), (5000, 5600), (5000, 5600), (5010, 5300)]); assert sorted(fp) == [((), 5000, 5600), (A, 0, 900)], fp  # inside an exon: dropped; the other group: one isoform
        import io; pr = io.StringIO(); c, n, a = assign(sam, pr); assert pr.getvalue() == 'r1\t2\nr2\t2\nr4\t1\n', pr.getvalue(); assert (c['2'], c['1'], n, a) == (2, 1, 4, 3), (c, n, a)  # r1 tie -> shorter transcript; r2 long deletion dropped on 1; r4 tie -> less clipping
        buf = io.StringIO(); gtf({'7': ('c', '+', 0, 900, A), '8': ('c', '.', 5000, 5600, ()), '9': ('c', '-', 0, 900, L)}, {'7': 3, '8': 4, '9': 2}, buf); g = buf.getvalue().splitlines()
        assert len(g) == 8 and g[0].split('\t')[:8] == ['c', 'flairlike', 'transcript', '1', '900', '.', '+', '.'] and g[1].split('\t')[3:5] == ['1', '100'] and g[5].split('\t')[3:5] == ['801', '900'] and g[6].split('\t')[2:5] == ['transcript', '5001', '5600'] and g[7].split('\t')[2] == 'exon', g  # 3-read floor drops 9
        print('ok'); sys.exit()
    if sys.argv[1] == 'firstpass':  # flair_like.py firstpass BAM GENOME PREFIX [THREADS] -> PREFIX.fa, PREFIX.tsv (id, contig, chain)
        import pysam
        bam, fa, pre = sys.argv[2:5]; g = pysam.FastaFile(fa); i = 0; e4 = int(__import__('os').environ.get('FLAIR_ENDS', '0'))  # 1 = E4, 2 = E5 (+ single-exon candidates)
        if e4: globals()['FP'] = e4  # set before the pool forks
        with multiprocessing.Pool(int(sys.argv[5]) if len(sys.argv) > 5 else 32) as p, open(pre + '.fa', 'w') as w, open(pre + '.tsv', 'w') as t:
            for contig, out in p.imap(work, [(bam, c) for c in pysam.AlignmentFile(bam).references]):  # imap: candidate ids and order fixed, so reruns agree
                for ch, a, b in (out if e4 else [(ch, v[4], v[5]) for ch, v in out.items() if not v[2]]):
                    i += 1; seq = ''.join(g.fetch(contig, x, y) for x, y in zip((a,) + tuple(k[1] for k in ch), tuple(k[0] for k in ch) + (b,)))
                    m = {{'GTAG': '+', 'CTAC': '-'}.get((g.fetch(contig, x, x + 2) + g.fetch(contig, y - 2, y)).upper()) for x, y in ch} - {None}
                    w.write(f'>{i}\n{seq}\n'); t.write(f'{i}\t{contig}\t' + ','.join(f'{x}-{y}' for x, y in ch) + f'\t{a}\t{b}\t{m.pop() if len(m) == 1 else "."}\n')
        sys.exit()
    if sys.argv[1] == 'place':  # minimap2 -a ... | flair_like.py place PREFIX.tsv OUT: read -> assigned isoform, as a molecules-like table
        ids = {}                # (r1645; strict_truth.py reads column 1 = read, 5 = contig, 16 = position)
        for l in open(sys.argv[2]): i, c, ch = l.rstrip('\n').split('\t')[:3]; ids[i] = (c, int(ch.split('-')[0]) if ch else 0)
        import io
        buf = io.StringIO(); cnt, reads, assigned = assign(sys.stdin, buf)
        with open(sys.argv[3], 'w') as w:
            for l in buf.getvalue().splitlines():
                q, i = l.split('\t'); c, pos = ids[i]; w.write('\t'.join([q] + ['NA'] * 3 + [c] + ['NA'] * 10 + [str(pos)]) + '\n')
        print(f'reads {reads}; assigned {assigned}'); sys.exit()
    if sys.argv[1] == 'gtf':  # minimap2 -a ... | flair_like.py gtf PREFIX.tsv OUT.gtf (driver stage flairlike, register 1640 E5)
        ids = {}
        for l in open(sys.argv[2]):
            i, c, ch, a, b, st = l.rstrip('\n').split('\t'); ids[i] = (c, st, int(a), int(b), tuple(tuple(map(int, x.split('-'))) for x in ch.split(',')) if ch else ())
        cnt, reads, assigned = assign(sys.stdin)
        with open(sys.argv[3], 'w') as w: gtf(ids, cnt, w)
        print(f'reads {reads}; assigned {assigned}; candidates {len(ids)}; isoforms {sum(n >= MIN for n in cnt.values())}'); sys.exit()
    if sys.argv[1] == 'assign':  # minimap2 -a ... | flair_like.py assign PREFIX.tsv REFSEQ DEFAULT_GTF FLAIR_GTF: FLAIR's read assignment, then the score
        ids = {}
        for l in open(sys.argv[2]):
            i, c, ch = l.rstrip('\n').split('\t')[:3]
            if ch: ids[i] = (c, tuple(tuple(map(int, x.split('-'))) for x in ch.split(',')))  # single-exon candidates are not scored
        cnt, reads, assigned = assign(sys.stdin)
        ref, dflt, fl = ({(c, ch) for c, s, ch in read(p)} for p in sys.argv[3:6])
        E3 = {ids[i] for i, n in cnt.items() if n >= MIN and i in ids}
        if len(sys.argv) > 6:  # per first-pass chain: id, contig, introns, assigned reads, in refseq, in flair
            with open(sys.argv[6], 'w') as w:
                for i, k in ids.items(): w.write(f'{i}\t{k[0]}\t{len(k[1])}\t{cnt.get(i, 0)}\t{int(k in ref)}\t{int(k in fl)}\n')
        print(f'reads {reads}; assigned {assigned}; first-pass chains {len(ids)}; FLAIR chains {len(fl)}\narm\tchains\tcorrect\tprecision\tshared\treached\tinside\tjaccard')
        for name, S in (('first pass', set(ids.values())), ('E3', E3), ('E3 + default', E3 | dflt)):
            n = len(S); ok = len(S & ref); sh = len(S & fl); print(f'{name}\t{n}\t{ok}\t{ok / max(n, 1):.3f}\t{sh}\t{sh / len(fl):.3f}\t{sh / max(n, 1):.3f}\t{sh / len(S | fl):.3f}')
        sys.exit()
    bam, refp, dp, flp, outp = sys.argv[1:6]; thr = int(sys.argv[6]) if len(sys.argv) > 6 else 32
    import pysam
    ref, dflt, fl = ({(c, ch) for c, s, ch in read(p)} for p in (refp, dp, flp)); arms = {'E0': set(), 'E1': set(), 'E2': set()}
    with multiprocessing.Pool(thr) as p, open(outp, 'w') as w:
        for contig, out in p.imap_unordered(work, [(bam, c) for c in pysam.AlignmentFile(bam).references]):
            for ch, (n, a1, dr, a2, _, _) in out.items():
                k = (contig, ch); w.write(f'{contig}\t{len(ch)}\t{n}\t{a1}\t{dr}\t{a2}\t{int(k in ref)}\t{int(k in fl)}\n')
                for name, x in (('E0', n), ('E1', a1), ('E2', a2)):
                    if x >= MIN: arms[name].add(k)
    arms['default'] = dflt; arms['E2 + default'] = arms['E2'] | dflt
    print(f'FLAIR chains {len(fl)}\narm\tchains\tcorrect\tprecision\tshared\treached\tinside\tjaccard')
    for name, S in arms.items():
        n = len(S); ok = len(S & ref); sh = len(S & fl); print(f'{name}\t{n}\t{ok}\t{ok / max(n, 1):.3f}\t{sh}\t{sh / len(fl):.3f}\t{sh / max(n, 1):.3f}\t{sh / len(S | fl):.3f}')
