#!/usr/bin/env python3
"""Cap signature of Iso-Seq reads: a non-templated G at the read's 5' end (register 1717).

  capg.py gate REF_GTF BAM FASTA [THREADS] [CONTIG,...]     share of capped reads at RefSeq 5' ends vs inside internal exons
  capg.py --selftest
A read is `cap` when its 5' soft clip is 1-3 bases, all G; `uninf` when the genome base just 5' of the alignment is G
(the extra G would align); else `no`. Read strand is taken as the transcript strand (FLNC reads)."""
import sys, collections, bisect, multiprocessing
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')

def cap_state(rev, cigar, seq, up):
    """rev: read on the minus strand; cigar: pysam tuples; up: genome base just 5' of the read in transcript sense."""
    op, n = cigar[-1] if rev else cigar[0]
    clip = (seq[len(seq) - n:].translate(COMP)[::-1] if rev else seq[:n]) if op == 4 else ''
    if 0 < len(clip) <= 3 and set(clip) == {'G'}: return 'cap'
    return 'uninf' if up.upper() == 'G' else 'no'

def read5(r, seq):
    """-> (strand, 5' position 0-based, state) of a pysam read on contig sequence seq."""
    if r.is_reverse:
        p = r.reference_end - 1; up = seq[p + 1:p + 2].translate(COMP)
    else:
        p = r.reference_start; up = seq[p - 1:p] if p else ''
    return '-' if r.is_reverse else '+', p, cap_state(r.is_reverse, r.cigartuples, r.query_sequence, up)

def gate_work(job):
    bam, fa, contig, tss, internal = job; import pysam
    seq = pysam.FastaFile(fa).fetch(contig); T = collections.Counter()
    near = lambda s, p, d: (lambda v, i: (i < len(v) and v[i] - p <= d) or (i and p - v[i - 1] <= d))(tss[s], bisect.bisect_left(tss[s], p))
    def inside(s, p):
        v = internal[s]; i = bisect.bisect_right(v, (p, 1 << 60)) - 1
        return i >= 0 and v[i][0] <= p < v[i][1]
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary: continue
        s, p, st = read5(r, seq)
        if 'N' in r.cigarstring and r.has_tag('ts'): T['ts', r.get_tag('ts')] += 1
        g = 'tss' if near(s, p, 10) else ('internal' if inside(s, p) and not near(s, p, 500) else None)
        if g: T[g, 'N' in r.cigarstring, st] += 1
    return T

def gate(ref, bam, fa, thr, only):
    import pysam, gzip, re
    ex = collections.defaultdict(list)
    for l in (gzip.open(ref, 'rt') if ref.endswith('.gz') else open(ref)):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
    tss = collections.defaultdict(lambda: {'+': set(), '-': set()}); inner = collections.defaultdict(lambda: {'+': [], '-': []})
    for (c, s, _), e in ex.items():
        e.sort(); tss[c][s].add(e[0][0] if s == '+' else e[-1][1] - 1); inner[c][s] += e[1:-1]
    have = set(pysam.AlignmentFile(bam).references)
    jobs = [(bam, fa, c, {s: sorted(v) for s, v in tss[c].items()}, {s: sorted(v) for s, v in inner[c].items()}) for c in tss if c in have and (not only or c in only)]
    T = collections.Counter()
    with multiprocessing.Pool(thr) as pool:
        for t in pool.imap_unordered(gate_work, jobs): T.update(t)
    print(f"# {bam}: spliced reads ts:+ {T['ts', '+']} ts:- {T['ts', '-']}")
    print('group\tspliced\treads\tinformative share\tcapped among informative')
    for g in ('tss', 'internal'):
        for sp in (True, False):
            c, n, u = (T[g, sp, k] for k in ('cap', 'no', 'uninf')); print(f"{g}\t{int(sp)}\t{c + n + u}\t{(c + n) / max(1, c + n + u):.3f}\t{c / max(1, c + n):.3f}")

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        assert cap_state(False, [(4, 1), (0, 5)], 'GACGTA', 'A') == 'cap' and cap_state(False, [(0, 6)], 'GACGTA', 'G') == 'uninf'
        assert cap_state(True, [(0, 5), (4, 2)], 'ACGTACC', 'T') == 'cap' and cap_state(True, [(0, 5), (4, 2)], 'ACGTACT', 'T') == 'no'
        assert cap_state(False, [(4, 4), (0, 2)], 'GGGGAC', 'G') == 'uninf'
        print('selftest ok')
    elif sys.argv[1] == 'gate':
        gate(*sys.argv[2:5], int(sys.argv[5]) if len(sys.argv) > 5 else 16, set(sys.argv[6].split(',')) if len(sys.argv) > 6 else None)

# ---- rule mode (PREREG row 1717) -------------------------------------------------------------------------------
def clusters(cap, inf):
    """cap / inf: position -> capped / informative reads. -> [(position, capped, informative)] of capped clusters, largest first."""
    c = dict(cap); out = []
    while c:
        m = max(c, key=lambda p: (c[p], -p)); win = [p for p in c if abs(p - m) <= 10]; n = sum(c[p] for p in win)
        ni = sum(v for p, v in inf.items() if abs(p - m) <= 10)
        for p in win: del c[p]
        if n >= 3 and n >= 0.5 * ni: out.append((m, n, ni))
    return sorted(out, key=lambda x: -x[1])

def rule_work(job):
    bam, fa_p, contig, D, ref, t5, peaks = job; import pysam
    from nested_sub import features
    from alt_ends import flags
    fa = pysam.FastaFile(fa_p); seq = fa.fetch(contig); dd = lambda: collections.defaultdict(collections.Counter)
    cnt = collections.Counter(); lo, hi, cap, inf = dd(), dd(), dd(), dd()
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for o, ln in r.cigartuples:
            if o == 3: ch.append((pos, pos + ln))
            if o in (0, 2, 3, 7, 8): pos += ln
        ch = tuple(ch); cnt[ch] += 1; lo[ch][r.reference_start] += 1; hi[ch][pos] += 1; _, p, st = read5(r, seq)
        if st != 'uninf': inf[ch][p] += 1; cap[ch][p] += st == 'cap'
    def pk(s, p): return int(peaks.dist(contig, s, p) == 0), int(peaks.dist(contig, s, p + (100 if s == '+' else -100)) == 0)
    best = {}
    for K, (strand, ks, ke) in D.items():
        for k in range(1, len(K)):
            for M, suf in ((K[k:], True), (K[:k], False)):
                if cnt.get(M, 0) >= 5 and cnt.get(K, 0) >= best.get((M, suf), (-1,))[0]: best[(M, suf)] = (cnt.get(K, 0), K, k, strand)
    rows = []
    for (M, suf), (rk, K, k, strand) in best.items():
        plus = strand == '+'; _, _, mode = features(M, K, k, suf, cnt, lo, hi); base = [contig, strand, int(M in D), int(M in ref.get(strand, ())), cnt[M], rk]
        if suf == plus:
            cl = clusters({p: v for p, v in cap[M].items() if v}, inf[M]); m5 = mode if plus else mode - 1
            rows.append(['NEST5'] + base + [int(M in t5), sum(inf[M].values()), len(cl), cl[0][1] if cl else 0, *(pk(strand, cl[0][0]) if cl else (0, 0)), *pk(strand, m5)])
        else:
            sig, arich = flags(fa.fetch, contig, plus, mode); rows.append(['NEST3'] + base + [int(arich and not sig)])
    for K, (strand, ks, ke) in D.items():
        if cnt.get(K, 0) < 5: continue
        plus = strand == '+'; c5 = lo[K] if plus else hi[K]; m = min(c5, key=lambda p: (-c5[p], p if plus else -p)); m5 = m if plus else m - 1
        cl = clusters({p: v for p, v in cap[K].items() if v}, inf[K]); alt = []
        for c in cl[1:]:
            if all(abs(c[0] - x[0]) >= 50 for x in [cl[0]] + alt): alt.append(c)
        ap = [pk(strand, c[0]) for c in alt]
        rows.append(['CHAIN', contig, strand, 1, int(K in ref.get(strand, ())), cnt[K], cnt[K], *pk(strand, ks if plus else ke - 1), *pk(strand, m5), len(cl),
                     *(pk(strand, cl[0][0]) if cl else (0, 0)), len(alt), sum(a for a, _ in ap), sum(b for _, b in ap)])
    return rows

def rule(ref_p, def_p, t5_p, bam, fa, bed, out, thr, skip):
    import pysam; from alt_events import load; from sq_nested import Peaks
    D = collections.defaultdict(dict); ref = collections.defaultdict(dict); t5 = collections.defaultdict(set); peaks = Peaks(bed)
    for (c, s), d in load(def_p).items():
        for ch, v in d.items(): D[c][ch] = (s, min(t[0] for t in v), max(t[1] for t in v))
    for (c, s), d in load(ref_p).items(): ref[c][s] = set(d)
    for (c, s), d in load(t5_p).items(): t5[c] |= set(d)
    have = set(pysam.AlignmentFile(bam).references)
    jobs = [(bam, fa, c, D[c], ref[c], t5[c], peaks) for c in sorted(D, key=lambda c: -len(D[c])) if c in have and (c in skip) == (thr < 0)]
    with multiprocessing.Pool(abs(thr)) as pool, open(out, 'w') as w:
        for rows in pool.imap_unordered(rule_work, jobs):
            for r in rows: w.write('\t'.join(map(str, r)) + '\n')
    score(out)

def score(path):
    R = collections.defaultdict(list)
    for l in open(path):
        f = l.rstrip('\n').split('\t'); R[f[0]].append([int(x) for x in f[3:]])
    sh = lambda a, b: f'{a} / {b} ({a / b:.3f})' if b else f'{a} / 0'
    # NEST5: emitted annotated reads k_reads trunc5 n_inf n_clusters top_cap clpk clbg mpk mbg
    N = R['NEST5']; ne = [r for r in N if not r[0]]; em = [r for r in N if r[0]]
    res = [r for r in ne if r[6]]; print(f'## A1 rescue: not emitted {len(ne)} (annotated {sum(r[1] for r in ne)}); with a capped cluster {len(res)}')
    print(f'   rescued annotated {sh(sum(r[1] for r in res), len(res))}; annotated-lost recovered {sh(sum(r[1] for r in res), sum(r[1] for r in ne))}')
    for name, S in (('rescued, cluster position', res), ('rescued unannotated', [r for r in res if not r[1]]), ('rescued annotated', [r for r in res if r[1]])):
        print(f'   {name}: in CAGE peak {sh(sum(r[8] for r in S), len(S))}; background {sh(sum(r[9] for r in S), len(S))}')
    S = [r for r in ne if not r[6]]; print(f'   not rescued, commonest start: in peak {sh(sum(r[10] for r in S), len(S))}; background {sh(sum(r[11] for r in S), len(S))}')
    ann = [r for r in em if r[1]]; un = [r for r in em if not r[1]]; drop = lambda r: r[5] >= 5 and not r[6]
    print(f'## A2 drop: emitted annotated {len(ann)}, unannotated {len(un)}')
    print(f'   capG rule: annotated kept {sh(sum(not drop(r) for r in ann), len(ann))}; unannotated dropped {sh(sum(drop(r) for r in un), len(un))}')
    print(f'   trunc5:    annotated kept {sh(sum(r[4] for r in ann), len(ann))}; unannotated dropped {sh(sum(not r[4] for r in un), len(un))}')
    for name, S in (('dropped unannotated', [r for r in un if drop(r)]), ('kept unannotated', [r for r in un if not drop(r)]), ('kept unannotated with a cluster', [r for r in un if r[6]]),
                    ('trunc5-dropped unannotated', [r for r in un if not r[4]]), ('dropped annotated', [r for r in ann if drop(r)]), ('kept annotated', [r for r in ann if not drop(r)])):
        print(f'   {name}, commonest start: in peak {sh(sum(r[10] for r in S), len(S))}; background {sh(sum(r[11] for r in S), len(S))}')
    # CHAIN: 1 annotated reads reads epk ebg mpk mbg n_clusters clpk clbg n_alt alt_pk alt_bg
    for name, S in (('annotated', [r for r in R['CHAIN'] if r[1]]), ('unannotated', [r for r in R['CHAIN'] if not r[1]])):
        C = [r for r in S if r[8]]; print(f'## B {name} emitted chains >= 5 reads: {len(S)}; with a capped cluster {sh(len(C), len(S))}')
        for lab, i in (('emitted 5\' end', 4), ('commonest start', 6), ('largest capped cluster', 9)):
            print(f'   B1 {lab} (chains with a cluster): in peak {sh(sum(r[i] for r in C), len(C))}; background {sh(sum(r[i + 1] for r in C), len(C))}')
        print(f'   B1 emitted 5\' end, chains without a cluster: in peak {sh(sum(r[4] for r in S if not r[8]), len(S) - len(C))}')
        n = sum(r[11] for r in S); print(f'   B2 further clusters >= 50 bp away: {n} on {sum(r[11] > 0 for r in S)} chains; in peak {sh(sum(r[12] for r in S), n)}; background {sh(sum(r[13] for r in S), n)}')
    T = [r for r in R['NEST3'] if r[0]]; ann = [r for r in T if r[1]]; dr = [r for r in T if r[4]]
    print(f'## C emitted early ends {len(T)}: annotated kept {sh(sum(not r[4] for r in ann), len(ann))}; dropped {len(dr)}, annotated among them {sh(sum(r[1] for r in dr), len(dr))}')

if __name__ == '__main__' and sys.argv[1] == 'rule':   # capg.py rule REF DEFAULT TRUNC5 BAM FASTA PEAKS OUT THREADS SKIP_CONTIGS (THREADS < 0: only those contigs)
    rule(*sys.argv[2:9], int(sys.argv[9]), set(sys.argv[10].split(',')))
elif __name__ == '__main__' and sys.argv[1] == 'score': score(sys.argv[2])

def ends(ref_p, a_p, b_p, bed, tsv, skip):   # PREREG row 1718: capg.py ends REF PREFIX.gtf PREFIX.final.gtf PEAKS RULE_TSV SKIP_CONTIGS
    from alt_events import load; from sq_nested import Peaks
    pk = Peaks(bed); ref = load(ref_p); A = load(a_p, 'reads'); B = load(b_p, 'reads'); T = collections.Counter()
    for (c, s), d in A.items():
        if c in skip: continue
        for ch, v in d.items():
            w = B.get((c, s), {}).get(ch)
            if not w or ch not in ref.get((c, s), ()) or max(t[2] or 0 for t in v) < 5: continue
            T['n'] += 1
            for name, m in (('most extended (PREFIX.gtf)', v), ('median (PREFIX.final.gtf)', w)):
                p = min(t[0] for t in m) if s == '+' else max(t[1] for t in m) - 1; q = p + (100 if s == '+' else -100)
                T[name, 'in'] += pk.dist(c, s, p) == 0; T[name, '50'] += pk.dist(c, s, p) <= 50; T[name, 'bg'] += pk.dist(c, s, q) == 0; T[name, 'bg50'] += pk.dist(c, s, q) <= 50
    print(f"annotated chains >= 5 reads in both files: {T['n']}")
    for name in ('most extended (PREFIX.gtf)', 'median (PREFIX.final.gtf)'):
        print(f"{name}: in peak {T[name, 'in'] / T['n']:.3f}; within 50 bp {T[name, '50'] / T['n']:.3f}; background {T[name, 'bg'] / T['n']:.3f} / {T[name, 'bg50'] / T['n']:.3f}")
    R = [l.split('\t') for l in open(tsv) if l.startswith('CHAIN')]; R = [r for r in R if r[4] == '1']
    print(f"commonest read start, row 1717 table, annotated chains >= 5 reads ({len(R)}): in peak {sum(int(r[9]) for r in R) / len(R):.3f}; emitted end there {sum(int(r[7]) for r in R) / len(R):.3f}")

if __name__ == '__main__' and sys.argv[1] == 'ends': ends(*sys.argv[2:7], set(sys.argv[7].split(',')))

# ---- stage mode (PREREG row 1720; driver stage `capstarts`) ---------------------------------------------------
def merged(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]: out[-1][1] = max(out[-1][1], b)
        else: out.append([a, b])
    return out

def stage_work(job):
    bam, fa_p, contig, D, ends = job; import pysam   # D: chain -> (strand, start, end, gene_id); ends: early-end mode (row 1722)
    from alt_ends import flags
    fa = pysam.FastaFile(fa_p); seq = fa.fetch(contig); end3 = collections.defaultdict(collections.Counter); dd = lambda: collections.defaultdict(collections.Counter)
    first = {'+': [], '-': []}; inner = {'+': [], '-': []}
    for K, (s, ks, ke, _) in D.items():
        first[s].append((ks, K[0][0]) if s == '+' else (K[-1][1], ke)); inner[s] += [(K[i][1], K[i + 1][0]) for i in range(len(K) - 1)]
    first = {s: merged(v) for s, v in first.items()}; inner = {s: merged(v) for s, v in inner.items()}
    def hit(v, p): i = bisect.bisect_right(v, [p, 1 << 60]) - 1; return i >= 0 and v[i][0] <= p < v[i][1]
    cnt = collections.Counter(); cap, inf = dd(), dd(); chk = collections.Counter()
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for o, ln in r.cigartuples:
            if o == 3: ch.append((pos, pos + ln))
            if o in (0, 2, 3, 7, 8): pos += ln
        ch = tuple(ch); cnt[ch] += 1; s, p, st = read5(r, seq)
        if ends: end3[ch][pos if s == '+' else r.reference_start] += 1; continue
        if st == 'uninf': continue
        inf[ch][p] += 1; cap[ch][p] += st == 'cap'; f, n = hit(first[s], p), hit(inner[s], p)
        if f != n: chk['first' if f else 'inner', st] += 1
    best = {}
    for K, (strand, ks, ke, gid) in D.items():
        for k in range(1, len(K)):
            M = K[k:] if (strand == '+') != ends else K[:k]
            if M not in D and cnt.get(M, 0) >= 5 and cnt.get(K, 0) >= best.get(M, (-1,))[0]: best[M] = (cnt.get(K, 0), strand, ks, ke, gid)
    lines = []
    for M, (_, strand, ks, ke, gid) in sorted(best.items()):
        if ends:   # commonest read end (ties to the more extended), peaked, polyA signal upstream, no A-rich run after
            c = end3[M]; plus = strand == '+'; m = min(c, key=lambda p: (-c[p], -p if plus else p)); sig, arich = flags(fa.fetch, contig, plus, m)
            if sum(v for p, v in c.items() if abs(p - m) <= 10) < 0.5 * cnt[M] or not sig or arich: continue
            cl = [(m, c[m])]; a, b = (ks, m) if plus else (m, ke)
        else:
            cl = clusters({p: v for p, v in cap[M].items() if v}, inf[M])
            if not cl: continue
            a, b = (cl[0][0], ke) if strand == '+' else (ks, cl[0][0] + 1)
        if not (a < M[0][0] and b > M[-1][1]): continue
        tid = f"{'END' if ends else 'CAP'}_{contig}_{cl[0][0]}_{len(M)}_{len(lines)}"; at = f'gene_id "{gid}"; transcript_id "{tid}";'
        ex = list(zip([a] + [i[1] for i in M], [i[0] for i in M] + [b]))
        lines.append(f'{contig}\trustle\ttranscript\t{a + 1}\t{b}\t.\t{strand}\t.\t{at} reads "{cnt[M]}"; ' + (f'end_reads "{cl[0][1]}"; earlyend "true";\n' if ends else f'cap_reads "{cl[0][1]}"; capstart "true";\n')
                     + ''.join(f'{contig}\trustle\texon\t{x + 1}\t{y}\t.\t{strand}\t.\t{at} exon_number "{n + 1}";\n' for n, (x, y) in enumerate(ex)))
    return chk, lines

def stage(gtf, bam, fa, out, thr, check_only, ends=False):
    import pysam, re
    ex = collections.defaultdict(list)
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1), re.search(r'gene_id "([^"]+)"', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
    D = collections.defaultdict(dict)
    for (c, s, _, g), e in ex.items():
        e.sort()
        if len(e) > 1:
            ch = tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)); o = D[c].get(ch)
            D[c][ch] = (s, min(e[0][0], o[1]) if o else e[0][0], max(e[-1][1], o[2]) if o else e[-1][1], g)
    have = set(pysam.AlignmentFile(bam).references); chk = collections.Counter(); new = []
    with multiprocessing.Pool(thr) as pool:
        for c, l in pool.imap(stage_work, [(bam, fa, c, D[c], ends) for c in sorted(D) if c in have]): chk.update(c); new += l
    sh = lambda g: chk[g, 'cap'] / max(1, chk[g, 'cap'] + chk[g, 'no']); ok = sh('first') >= 0.5 and sh('inner') <= 0.15
    msg = f"self-check: capped share of reads starting in a first exon {sh('first'):.3f} (needs >= 0.5), inside an internal exon {sh('inner'):.3f} (needs <= 0.15)"
    if ends: msg = 'early ends: peaked read end, polyA signal, no A-rich run'
    elif not ok: sys.exit(f'capstarts refused: this library has no usable cap G; {msg}')
    if not check_only:
        with open(out, 'w') as w:
            for l in open(gtf): w.write(l)
            w.writelines(new)
    print(f"{msg}; added {len(new)} nested {'early ends' if ends else 'starts with a capped cluster'}")

def exempt_work(job):
    bam, fa_p, contig, C, sub, other, known = job; import pysam   # C: candidate chain -> (strand, annotated); sub: terminal sub-runs of default chains; other: second BAM or None; known: annotated introns
    from alt_ends import flags
    fa = pysam.FastaFile(fa_p); seq = fa.fetch(contig); dd = lambda: collections.defaultdict(collections.Counter)
    cnt = collections.Counter(); cap, inf, end3, st5 = dd(), dd(), dd(), dd()
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for o, ln in r.cigartuples:
            if o == 3: ch.append((pos, pos + ln))
            if o in (0, 2, 3, 7, 8): pos += ln
        ch = tuple(ch)
        if ch not in C: continue
        cnt[ch] += 1; s, p, st = read5(r, seq); end3[ch][pos if s == '+' else r.reference_start] += 1; st5[ch][p] += 1
        if st != 'uninf': inf[ch][p] += 1; cap[ch][p] += st == 'cap'
    oc = collections.Counter()   # row 1733: full-chain primary reads of the candidate in an independent sample
    for r in (pysam.AlignmentFile(other).fetch(contig) if other else ()):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for o, ln in r.cigartuples:
            if o == 3: ch.append((pos, pos + ln))
            if o in (0, 2, 3, 7, 8): pos += ln
        if tuple(ch) in C: oc[tuple(ch)] += 1
    rows = []
    for M, (strand, ann) in C.items():
        pa = 0; m = -1
        if cnt[M]:   # the earlyends 3' rule: commonest read end, peaked, polyA signal upstream, no A-rich run after
            c = end3[M]; plus = strand == '+'; m = min(c, key=lambda p: (-c[p], -p if plus else p)); sig, arich = flags(fa.fetch, contig, plus, m)
            pa = int(sum(v for p, v in c.items() if abs(p - m) <= 10) >= 0.5 * cnt[M] and sig and not arich)
        cl = clusters({p: v for p, v in cap[M].items() if v}, inf[M])   # last two columns (row 1731): 5' position (capped cluster, else commonest read start), commonest 3' end
        rows.append((contig, strand, len(M), cnt[M], sum(inf[M].values()), int(bool(cl)), pa, int(M in sub), ann, cl[0][0] if cl else (st5[M].most_common(1)[0][0] if cnt[M] else -1), m, int(all(i in known for i in M)), oc[M] if other else -1))
    return rows

def exempt(ref_p, def_p, off_p, bam, fa, out, thr, other=None):   # PREREG row 1724: capg.py exempt REF DEFAULT_GTF FRACTION_OFF_GTF BAM FASTA OUT_TSV THREADS
    import pysam; from alt_events import load
    D, N, R = load(def_p), load(off_p), load(ref_p); have = set(pysam.AlignmentFile(bam).references); jobs = collections.defaultdict(lambda: ({}, set(), set()))
    tot = collections.Counter()
    for k in set(D) | set(N):
        d, n, r = set(D.get(k, ())), set(N.get(k, ())), set(R.get(k, ()))
        tot['default'] += len(d); tot['default_ok'] += len(d & r); tot['gone'] += len(d - n); tot['gone_ok'] += len((d - n) & r)
        if k[0] not in have: continue
        C, sub, known = jobs[k[0]]; known.update(i for ch in r for i in ch)
        for ch in n - d: C[ch] = (k[1], int(ch in r))
        for K in d: sub.update(K[i:] for i in range(1, len(K))); sub.update(K[:i] for i in range(1, len(K)))
    with multiprocessing.Pool(thr) as pool, open(out, 'w') as w:
        w.write(f"#default\t{tot['default']}\t{tot['default_ok']}\n#gone\t{tot['gone']}\t{tot['gone_ok']}\n")
        for rows in pool.imap_unordered(exempt_work, [(bam, fa, c, C, sub, other, known) for c, (C, sub, known) in sorted(jobs.items(), key=lambda x: -len(x[1][0]))]):
            for r in rows: w.write('\t'.join(map(str, r)) + '\n')
    print(exempt_score(open(out)))

def exempt_score(lines):
    """Rows: contig strand introns reads inf5 cap5 pa3 sub annotated. Decision set E = cap5 & pa3 & not a terminal sub-run."""
    H = {}; X = []
    for l in lines:
        f = l.split('\t')
        if l[0] == '#': H[f[0]] = (int(f[1]), int(f[2]))
        else: X.append(tuple(map(int, f[2:])))
    sh = lambda v: f'{len(v)}\t{sum(x[6] for x in v)}\t{sum(x[6] for x in v) / max(1, len(v)):.3f}'
    d, dk = H['#default']; P = dk / d; out = [f'default {d} chains, {dk} annotated ({P:.3f}); default chains gone with the rule off: {H["#gone"][0]} / annotated {H["#gone"][1]}',
                                               f'candidates (fraction-off, not in default): {sh(X)}', 'group\tn\tannotated\tshare']
    und = [x for x in X if x[2] == 0]; dec = [x for x in X if x[2] > 0]; both = [x for x in dec if x[3] and x[4]]; E = [x for x in both if not x[5]]
    rest = sorted((x for x in dec if not (x[3] and x[4]) and not x[5]), key=lambda x: -x[1]); ctl = rest[:len(E)]
    for name, v in (('both ends', both), ('both ends, terminal sub-run', [x for x in both if x[5]]), ('E: both ends, not a sub-run', E),
                    ('control: top-read no-signature, not a sub-run, same n', ctl), ('cap only', [x for x in dec if x[3] and not x[4]]),
                    ('polyA end only', [x for x in dec if x[4] and not x[3]]), ('neither', [x for x in dec if not x[3] and not x[4]]),
                    ('undecidable 5\' (no informative read)', und), ('secondary: polyA end, any 5\', not a sub-run', [x for x in X if x[4] and not x[5]])):
        out.append(f'{name}\t{sh(v)}')
    ek = sum(x[6] for x in E); es = ek / max(1, len(E)); cs = sum(x[6] for x in ctl) / max(1, len(ctl)); comb = (dk + ek) / (d + len(E))
    for name, ok, txt in (('B1 quality', es >= P, f'{es:.3f} vs default {P:.3f}'), ('B2 depth control', es >= 1.5 * cs, f'{es:.3f} vs 1.5 x {cs:.3f}'),
                          ('B3 yield', ek >= 0.01 * dk, f'{ek} vs {0.01 * dk:.0f}'), ('B4 combined precision', comb >= P - 0.005, f'{comb:.4f} vs {P - 0.005:.4f}')):
        out.append(f"{name}\t{'PASS' if ok else 'FAIL'}\t{txt}")
    return '\n'.join(out)

def exempt_junc(rows, seed=0):
    """PREREG row 1733. Chain-level yardsticks of the unannotated decision set: all introns annotated; chain spelled by a read of the other sample."""
    import random
    X = [tuple(map(int, f[2:])) for f in rows if f[0][0] != '#']; X = [x for x in X if x[1] > 0 and not x[5]]   # with reads, not a sub-run
    sig = lambda x: x[2] > 0 and x[3] and x[4]; U = [x for x in X if sig(x) and not x[6]]; A = [x for x in X if sig(x) and x[6]]
    pool = collections.defaultdict(list); need = collections.Counter(x[1] for x in U); ctl = []
    for x in X:
        if x[2] > 0 and not sig(x) and not x[6]: pool[x[1]].append(x)
    rng = random.Random(seed)
    for n, k in sorted(need.items()): ctl += rng.sample(pool[n], min(k, len(pool[n])))
    G = {'E unannotated': U, 'E annotated': A, 'control unannotated (no signature, read-count matched)': ctl}
    r = {k: (sum(x[9] for x in v) / max(1, len(v)), sum(x[10] > 0 for x in v) / max(1, len(v)), sum(x[10] > 0 for x in v if not x[9]) / max(1, sum(1 for x in v if not x[9]))) for k, v in G.items()}
    out = ['group\tn\tall introns annotated\tchain read in the other sample\tsame, chains with a novel intron'] + [f'{k}\t{len(G[k])}\t{r[k][0]:.3f}\t{r[k][1]:.3f}\t{r[k][2]:.3f}' for k in G]
    u, a, c = (r[k] for k in G)
    for name, ok, txt in (('B1 known junctions vs control', u[0] >= 1.5 * c[0], f'{u[0]:.3f} vs 1.5 x {c[0]:.3f}'), ('B2 replication vs annotated', u[1] >= 0.5 * a[1], f'{u[1]:.3f} vs 0.5 x {a[1]:.3f}'),
                          ('B3 replication vs control', u[1] >= 1.5 * c[1], f'{u[1]:.3f} vs 1.5 x {c[1]:.3f}')):
        out.append(f"{name}\t{'PASS' if ok else 'FAIL'}\t{txt}")
    return '\n'.join(out)

if __name__ == '__main__' and sys.argv[1] == 'exempt': exempt(*sys.argv[2:8], int(sys.argv[8]), *sys.argv[9:10]); sys.exit(0)   # optional 9th argument: OTHER_BAM
if __name__ == '__main__' and sys.argv[1] == 'exempt-junc': print(exempt_junc(l.rstrip('\n').split('\t') for l in open(sys.argv[2]))); sys.exit(0)
def exempt_ends(rows, cage, near):
    """PREREG row 1731. rows: exempt TSV fields; cage(contig, strand, p) / near(contig, strand, p) -> supported. Text table + bars."""
    X = [(f[0], f[1]) + tuple(map(int, f[2:])) for f in rows if f[0][0] != '#']; X = [x for x in X if x[9] >= 0 and not x[7]]   # with reads, not a sub-run
    sig = lambda x: x[4] > 0 and x[5] and x[6]; E = [x for x in X if sig(x)]; n = sum(1 for x in E if not x[8])
    ctl = sorted((x for x in X if x[4] > 0 and not sig(x) and not x[8]), key=lambda x: -x[3])[:n]
    def rate(v):
        r = collections.Counter()
        for x in v:
            c, s, p5, p3 = x[0], x[1], x[9], x[10]
            d = 1 if s == '+' else -1; a, b = cage(c, s, p5), near(c, s, p3)
            r['5'] += a; r['3'] += b; r['both'] += a and b; r['5bg'] += cage(c, s, p5 + 100 * d); r['3bg'] += near(c, s, p3 - 200 * d)
        return {k: r[k] / max(1, len(v)) for k in ('5', '5bg', '3', '3bg', 'both')}
    G = {'E unannotated': [x for x in E if not x[8]], 'E annotated': [x for x in E if x[8]], 'control unannotated (top-read, no signature, same n)': ctl}
    R = {k: rate(v) for k, v in G.items()}; out = ['group\tn\t5\' in CAGE peak\tbackground\t3\' at polyA site\tbackground\tboth']
    out += [f"{k}\t{len(G[k])}\t" + '\t'.join(f'{R[k][c]:.3f}' for c in ('5', '5bg', '3', '3bg', 'both')) for k in G]
    u, a, c = R['E unannotated'], R['E annotated'], R['control unannotated (top-read, no signature, same n)']
    for name, ok, txt in (("B1 5' end", u['5'] >= 0.5 * a['5'], f"{u['5']:.3f} vs 0.5 x {a['5']:.3f}"), ("B2 3' end", u['3'] >= 0.5 * a['3'], f"{u['3']:.3f} vs 0.5 x {a['3']:.3f}"),
                          ('B3 both ends vs control', u['both'] >= 2 * c['both'], f"{u['both']:.3f} vs 2 x {c['both']:.3f}")):
        out.append(f"{name}\t{'PASS' if ok else 'FAIL'}\t{txt}")
    return '\n'.join(out)

if __name__ == '__main__' and sys.argv[1] == 'exempt-ends':   # capg.py exempt-ends EXEMPT_TSV CAGE_BED POLYA_BED
    from sq_nested import Peaks; pk = Peaks(sys.argv[3]); S = collections.defaultdict(list)
    for l in open(sys.argv[4]): f = l.split('\t'); S[f[0], f[5].strip()].append(int(f[1]))
    for v in S.values(): v.sort()
    def near(c, s, p): v = S.get((c, s), ()); i = bisect.bisect_left(v, p - 25); return int(i < len(v) and v[i] <= p + 25)
    print(exempt_ends((l.rstrip('\n').split('\t') for l in open(sys.argv[2])), lambda c, s, p: int(pk.dist(c, s, p) == 0), near)); sys.exit(0)
if __name__ == '__main__' and sys.argv[1] == 'exempt-selftest':
    r3 = lambda reads, sg, ann, allann, other: ['c', '+', '3', str(reads), '9', str(sg), str(sg), '0', str(ann), '5', '9', str(allann), str(other)]
    t = exempt_junc([r3(7, 1, 0, 1, 2), r3(7, 1, 0, 0, 0), r3(7, 1, 1, 1, 4), r3(7, 0, 0, 0, 0), r3(7, 0, 0, 0, 3), r3(90, 0, 0, 1, 9)])
    assert 'E unannotated\t2\t0.500\t0.500\t0.000' in t and 'read-count matched)\t2\t0.000\t0.500' in t and 'B1 known junctions vs control\tPASS' in t and 'B3 replication vs control\tFAIL' in t, t
    r2 = lambda reads, cap5, pa3, ann, p5, p3: ['c', '+', '3', str(reads), '9', str(cap5), str(pa3), '0', str(ann), str(p5), str(p3)]
    t = exempt_ends([r2(9, 1, 1, 0, 1000, 5000), r2(9, 1, 1, 0, 2000, 6000), r2(9, 1, 1, 1, 1000, 5000), r2(50, 0, 0, 0, 3000, 7000), r2(40, 0, 1, 0, 1000, 5000), r2(1, 0, 0, 0, 1000, 5000)],
                    lambda c, s, p: int(p == 1000), lambda c, s, p: int(p == 5000))
    assert 'E unannotated\t2\t0.500\t0.000\t0.500\t0.000\t0.500' in t and 'same n)\t2\t0.500\t0.000\t0.500\t0.000\t0.500' in t and 'B3 both ends vs control\tFAIL' in t and "B1 5' end\tPASS" in t, t
    row = lambda reads, inf5, cap5, pa3, sub, ann: f'c\t+\t3\t{reads}\t{inf5}\t{cap5}\t{pa3}\t{sub}\t{ann}\n'
    t = exempt_score(['#default\t1000\t200\n', '#gone\t5\t1\n'] + [row(9, 9, 1, 1, 0, 1)] * 3 + [row(9, 9, 1, 1, 0, 0)] + [row(9, 9, 1, 1, 1, 1)]
                     + [row(50, 50, 0, 0, 0, 0)] * 6 + [row(4, 0, 0, 1, 0, 1)])
    assert 'E: both ends, not a sub-run\t4\t3\t0.750' in t and 'control: top-read no-signature, not a sub-run, same n\t4\t0\t0.000' in t, t
    assert "undecidable 5' (no informative read)\t1\t1" in t and 'B1 quality\tPASS' in t and 'B3 yield\tPASS\t3 vs 2' in t and 'both ends, terminal sub-run\t1\t1' in t, t
    print('exempt selftest ok'); sys.exit(0)

if __name__ == '__main__' and sys.argv[1] == 'stage':   # capg.py stage GTF BAM FASTA OUT_GTF [THREADS] [--check-only] [--ends]
    stage(*sys.argv[2:6], int(sys.argv[6]) if len(sys.argv) > 6 else 16, '--check-only' in sys.argv, '--ends' in sys.argv)
