#!/usr/bin/env python3
"""PREREG r1637 (docs/PREREG_subset_merge_2026-11-11.md): merge reads that are subsets of one longer read.

  subset_merge.py BAM REFSEQ DEFAULT_GTF STRINGTIE_GTF FLAIR_GTF OUT_TSV [THREADS]      (selftest: subset_merge.py selftest)"""
import sys, collections, multiprocessing
sys.path.insert(0, __import__('os').path.dirname(__file__))
from tool_match import read
TOL = 10

def merge(reads):
    """reads: [(start, end, chain)] of one contig. Returns {chain: [full, merged_unique, merged_any, end_shared]}."""
    info = {}  # chain -> [count, min start, max end]
    for a, b, ch in reads:
        v = info.setdefault(ch, [0, a, b]); v[0] += 1; v[1] = min(v[1], a); v[2] = max(v[2], b)
    by = collections.defaultdict(set)
    for ch in info:
        for i in ch: by[i].add(ch)
    def contains(c, s, a, b):
        n = len(s)
        for k in range(len(c) - n + 1):
            if c[k:k + n] == s:
                lo = c[k - 1][1] if k else info[c][1]; hi = c[k + n][0] if k + n < len(c) else info[c][2]
                return a >= lo - TOL and b <= hi + TOL
        return False
    out = {ch: [v[0], 0, 0, 0] for ch, v in info.items()}
    for a, b, s in reads:
        cs = [c for c in by[s[0]] if len(c) > len(s) and contains(c, s, a, b)]
        if not cs: continue
        mx = [c for c in cs if not any(len(d) > len(c) and d != c and contains(d, c, info[c][1], info[c][2]) for d in cs)]
        for c in mx: out[c][2] += 1
        if len(mx) == 1:
            c = mx[0]; out[c][1] += 1; out[c][3] += abs(a - info[c][1]) <= TOL or abs(b - info[c][2]) <= TOL
    return out

def work(job):
    import pysam
    bam, contig = job; reads = []
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for op, ln in r.cigartuples:
            if op == 3: ch.append((pos, pos + ln))
            if op in (0, 2, 3, 7, 8): pos += ln
        reads.append((r.reference_start, pos, tuple(ch)))
    return contig, merge(reads)

if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        C = ((10, 20), (30, 40), (50, 60)); D = ((10, 20), (30, 40), (50, 60), (70, 80)); E = ((10, 20), (30, 45))
        r = merge([(0, 70, C), (22, 48, ((30, 40),)), (0, 28, ((10, 20),)), (5, 48, ((30, 40),)), (0, 100, E)])
        assert r[C] == [1, 1, 2, 0], r[C]   # (22,48) merged; (0,28,(10,20)) fits C and E: two maximal containers; (5,..) starts before the exon
        r = merge([(0, 70, C), (0, 100, D), (22, 48, ((30, 40),))]); assert r[D][1] == 2 and r[C][1] == 0, r  # nested containers: the short read and the C read both go to the outer one
        print('ok'); sys.exit()
    bam, refp, dp, stp, flp, outp = sys.argv[1:7]; thr = int(sys.argv[7]) if len(sys.argv) > 7 else 32
    import pysam
    ref, dflt, st, fl = (read(p) for p in (refp, dp, stp, flp)); strip = lambda S: {(c, ch) for c, s, ch in S}
    ref, dflt, st, fl = map(strip, (ref, dflt, st, fl))
    contigs = [c for c in pysam.AlignmentFile(bam).references]
    rows = collections.defaultdict(lambda: [0, 0, 0, 0]); merged = shared = 0
    with multiprocessing.Pool(thr) as p, open(outp, 'w') as w:
        for contig, out in p.imap_unordered(work, [(bam, c) for c in contigs]):
            for ch, (full, mu, ma, es) in out.items():
                k = (contig, ch); merged += mu; shared += es
                # every read chain (r1638): contig, introns, full, merged unique, merged any, in default/refseq/stringtie/flair
                w.write(f'{contig}\t{len(ch)}\t{full}\t{mu}\t{ma}\t{int(k in dflt)}\t{int(k in ref)}\t{int(k in st)}\t{int(k in fl)}\n')
                if k in dflt or full != 1: continue
                for name, m in (('unique container', mu), ('any maximal container', ma)):
                    if m >= 1: v = rows[name]; v[0] += 1; v[1] += k in ref; v[2] += k in st; v[3] += k in fl
    print(f'sets\tdefault {len(dflt)}\trefseq {len(ref)}\tstringtie {len(st)}\tflair {len(fl)}')
    print('rule\tnew_chains\tcorrect\tshare\tin_stringtie\tin_flair')
    for name, (n, ok, a, b) in rows.items(): print(f'{name}\t{n}\t{ok}\t{ok / max(n, 1):.3f}\t{a}\t{b}')
    print(f'merged reads {merged}; ending within {TOL} bp of an outer end of their container {shared} = {shared / max(merged, 1):.3f}')
