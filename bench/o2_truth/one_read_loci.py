#!/usr/bin/env python3
"""One read as the seed of a locus: what the lone reads at empty places are (register 1987).

    one_read_loci.py CAND.gtf BAM AS_TABLE TRUTH.gtf|- [--replicates A.bam,B.bam] [--stringtie ST.gtf] [--base BASE.gtf]

CAND.gtf = tools/multicopy_rescue.py candidates BASE FLOOR1 (spliced transcripts of the one-read-floor assembly where the
main assembly has none on that strand). A candidate is taken when exactly ONE primary read has its intron chain (the lone
read). Variants: L0 every such candidate; L1 the lone read has no other alignment within 10 ms of its best (AS_TABLE);
L2 = L1 and >= 2 introns. Correct = exact intron chain of TRUTH.gtf. --replicates: a primary read with exactly the chain
in another library. --stringtie: StringTie's transcripts with the same chain. --base: the main assembly, for the totals."""
import sys, os, re, collections, multiprocessing
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import thin_isoforms as T
G = {}

def scan(contig):
    """lone-read lookup: {chain: [read names]} for the wanted chains on a contig (primary alignments)"""
    import pysam
    want = G['by'].get(contig)
    if not want: return {}
    out = collections.defaultdict(list)
    for r in pysam.AlignmentFile(G['bam']).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        c = []; p = r.reference_start
        for o, k in r.cigartuples:
            if o == 3: c.append((p, p + k))
            if o in (0, 2, 3, 7, 8): p += k
        c = tuple(c)
        if c in want and len(out[c]) < 3: out[c].append((r.query_name, r.mapping_quality))
    return {(contig, c): v for c, v in out.items()}

def reads_with(bam, chains, thr):
    import pysam
    by = collections.defaultdict(set)
    for c, s, ch in chains: by[c].add(ch)
    G.update(bam=bam, by=by)
    with multiprocessing.get_context('fork').Pool(thr) as pool: res = pool.map(scan, [c for c in pysam.AlignmentFile(bam).references if c in by], chunksize=1)
    out = {}
    for d in res: out.update(d)
    return out

def main(cand, bam, tab, truth, *a):
    opt = dict(zip(a[::2], a[1::2])); thr = int(opt.get('--threads', 32))
    C = T.chains_gtf(cand); R = reads_with(bam, C, thr)
    lone = {c: R[(c[0], c[2])][0] for c in C if len(R.get((c[0], c[2]), ())) == 1}
    names = {v[0] for v in lone.values()}; tie = {}
    for l in open(tab):
        if l[0] == '#': continue
        f = l.split('\t', 9)
        if f[0] in names: tie[f[0]] = (int(f[7]) >= 0 and int(f[7]) >= int(f[6]) - 10)
    tr = T.chains_gtf(truth) if truth != '-' else None
    rep = set()
    for b in filter(None, opt.get('--replicates', '').split(',')):
        rep |= {(k[0], k[1]) for k in reads_with(b, set(lone), thr)}
    st = T.chains_gtf(opt['--stringtie']) if '--stringtie' in opt else None
    V = {'L0 every lone read': lambda c: True, 'L1 placed without a rival': lambda c: not tie.get(lone[c][0], False), 'L2 and >= 2 introns': lambda c: not tie.get(lone[c][0], False) and len(c[2]) >= 2,
         'left out by L1 (the read has a rival within 10 ms)': lambda c: tie.get(lone[c][0], False), 'L1 with exactly 1 intron': lambda c: not tie.get(lone[c][0], False) and len(c[2]) == 1}
    print(f'candidates at empty places {len(C)}; with exactly one read of that chain {len(lone)}')
    print('variant\tloci added\tcorrect / annotated\tshare\treplicated in another library\tshare\tStringTie has the chain\tof those, correct')
    for lab, f in V.items():
        s = [c for c in lone if f(c)]; n = len(s); ok = sum(c in tr for c in s) if tr is not None else None; rp = sum((c[0], c[2]) in rep for c in s) if '--replicates' in opt else None
        sc = [c for c in s if st is not None and c in st]
        print(lab, n, ok if ok is not None else 'NA', f'{ok / max(1, n):.3f}' if ok is not None else 'NA', rp if rp is not None else 'NA', f'{rp / max(1, n):.3f}' if rp is not None else 'NA',
              len(sc) if st is not None else 'NA', f'{sum(c in tr for c in sc)} ({sum(c in tr for c in sc) / max(1, len(sc)):.3f})' if st is not None and tr is not None else 'NA', sep='\t')
    if '--base' in opt:
        B = T.chains_gtf(opt['--base']); line = f'main assembly: {len(B)} spliced chains'
        if tr is not None: line += f', {len(B & tr)} correct ({len(B & tr) / max(1, len(B)):.4f})'
        if '--replicates' in opt:
            rb = set()
            for b in opt['--replicates'].split(','): rb |= {(k[0], k[1]) for k in reads_with(b, B, thr)}
            line += f"; replicated in another library {sum((c[0], c[2]) in rb for c in B) / max(1, len(B)):.3f}"
        print(line)
        if st is not None and tr is not None:
            so = st - B; print(f'StringTie: {len(st)} spliced chains, {len(st & tr)} correct ({len(st & tr) / max(1, len(st)):.4f}); chains the main assembly lacks {len(so)}, correct {len(so & tr)} ({len(so & tr) / max(1, len(so)):.3f})')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: print('one_read_loci selftest ok (pure parts are thin_isoforms)')
    else: main(*sys.argv[1:])
