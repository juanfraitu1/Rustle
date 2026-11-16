#!/usr/bin/env python3
"""Addendum 4 of PREREG_ptr_ppa_frozen_rules_2026-11-15: do 5' shortened chains start at one position, in both species.

  xspecies_starts.py Q.gtf Q.refseq.gtf Q.bam MAP.paf LIFT.sh T.bam OUT.tsv   |   --selftest"""
import sys, re, collections, subprocess, pysam
W = 10

def models(p):
    ex = collections.defaultdict(list); sd = {}
    for l in open(p):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': k = (f[0], re.search(r'transcript_id "([^"]+)"', f[8]).group(1)); ex[k].append((int(f[3]), int(f[4]))); sd[k] = f[6]
    out = {}
    for k, e in ex.items():
        if len(e) > 1: e.sort(); out[(k[0], tuple((a[1], b[0] - 1) for a, b in zip(e, e[1:])))] = sd[k]
    return out

def klass(ch, s, refs):
    """Y exact; S5 / S3 = contiguous sub-chain of a reference chain missing introns at the 5' / 3' side only; else None."""
    best = None
    for b in refs:
        n = len(ch)
        for i in range(len(b) - n + 1):
            if b[i:i + n] == ch:
                if len(b) == n: return 'Y'
                left = i > 0; right = i + n < len(b); m5, m3 = (left, right) if s == '+' else (right, left)
                k = 'S5' if m5 and not m3 else 'S3' if m3 and not m5 else None
                if k and best is None: best = k
                elif k and best != k: best = 'both'   # ambiguous between references
    return best if best in ('S5', 'S3') else None

def conc(E, at=None):
    """share of ends within W of `at` (default: the commonest end), that position, n."""
    n = sum(E.values())
    if not n: return 0.0, None, 0
    m = max(E, key=lambda p: (E[p], -p)) if at is None else at
    return sum(c for p, c in E.items() if abs(p - m) <= W) / n, m, n

def ends(bam, want):
    """{(contig, chain): {'L': Counter(reference_start), 'R': Counter(reference_end)}} of primary reads with exactly a wanted chain."""
    out = collections.defaultdict(lambda: (collections.Counter(), collections.Counter()))
    for r in pysam.AlignmentFile(bam, threads=8):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        pos = r.reference_start; ch = []
        for op, n in r.cigartuples:
            if op == 3: ch.append((pos, pos + n))
            if op in (0, 2, 3, 7, 8): pos += n
        k = (r.reference_name, tuple(ch))
        if ch and k in want: out[k][0][r.reference_start] += 1; out[k][1][pos] += 1
    return out

def lift(paf, sh, bed_rows, path):
    with open(path, 'w') as o:
        for c, a, b in sorted(set(bed_rows)): o.write(f'{c}\t{a}\t{b}\n')
    L = {}
    for l in subprocess.run(['bash', sh, paf, path], capture_output=True, text=True, check=True).stdout.splitlines():
        f = l.split('\t'); c, a, b = f[3].rsplit('_', 2); L[(c, int(a), int(b))] = (f[0], int(f[1]), int(f[2]))
    return L

def main(qg, qr, qbam, paf, sh, tbam, out):
    Q = models(qg); R = models(qr); by = collections.defaultdict(list)
    for (c, ch), s in R.items(): by[(c, s, ch[0])].append(ch); [by[(c, s, i)].append(ch) for i in ch[1:]]
    cls = {}
    for (c, ch), s in Q.items():
        k = klass(ch, s, set(by[(c, s, ch[0])]))
        if k: cls[(c, ch)] = k
    E = ends(qbam, cls); own = {}
    for k, (Lc, Rc) in E.items():
        cc, m, n = conc(Lc if Q[k] == '+' else Rc)
        if n >= 6: own[k] = (cc, m, n)
    L = lift(paf, sh, [(c, a, b) for (c, ch) in own for a, b in ch] + [(c, own[(c, ch)][1], own[(c, ch)][1] + 1) for (c, ch) in own], out + '.bed')
    tgt = {}
    for (c, ch), (cc, m, n) in own.items():
        t = [L.get((c, a, b)) for a, b in ch]; p = L.get((c, m, m + 1))
        if None in t or p is None or len({x[0] for x in t} | {p[0]}) != 1: continue
        first = t[0] if Q[(c, ch)] == '+' else t[-1]; srt = tuple(sorted((a, b) for _, a, b in t))
        tgt[(c, ch)] = ((t[0][0], srt), 0 if (first[1], first[2]) == srt[0] else 1, p[1])   # which read end is the 5' end in the target
    TE = ends(tbam, {v[0] for v in tgt.values()}); res = collections.defaultdict(list)
    with open(out, 'w') as o:
        o.write('chrom\tclass\tn_introns\town_reads\town_conc\tsister_reads\tsister_conc\n')
        for k, (tk, side, p) in tgt.items():
            sc, _, sn = conc(TE[tk][side], p) if tk in TE else (0.0, None, 0)
            o.write(f'{k[0]}\t{cls[k]}\t{len(k[1])}\t{own[k][2]}\t{own[k][0]:.3f}\t{sn}\t{sc:.3f}\n')
            if sn >= 3: res[cls[k]].append((own[k][0], sc))
    n_all = collections.Counter(cls.values()); n_own = collections.Counter(cls[k] for k in own); y = None
    print('class\tchains\t>= 6 own reads\tqualify (>= 3 sister reads)\town conc >= 0.5\tsister conc >= 0.5\tstart-like (both)\tx Y')
    for k in ('Y', 'S3', 'S5'):
        v = res[k]; n = max(1, len(v)); sl = sum(a >= 0.5 and b >= 0.5 for a, b in v) / n; y = sl if k == 'Y' else y
        print(f"{k}\t{n_all[k]}\t{n_own[k]}\t{len(v)}\t{sum(a >= 0.5 for a, _ in v) / n:.3f}\t{sum(b >= 0.5 for _, b in v) / n:.3f}\t{sl:.3f}\t{sl / y if y else 0:.2f}")

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        r = {((1, 2), (3, 4), (5, 6))}
        assert klass(((1, 2), (3, 4), (5, 6)), '+', r) == 'Y' and klass(((3, 4), (5, 6)), '+', r) == 'S5' and klass(((3, 4), (5, 6)), '-', r) == 'S3'
        assert klass(((1, 2), (3, 4)), '+', r) == 'S3' and klass(((3, 4),), '+', r) is None and klass(((1, 2), (5, 6)), '+', r) is None
        assert klass(((3, 4), (5, 6)), '+', r | {((3, 4), (5, 6), (7, 8))}) is None      # 5' short of one reference, 3' short of another
        assert conc(collections.Counter({100: 5, 105: 2, 300: 3})) == (0.7, 100, 10) and conc(collections.Counter({100: 5, 300: 5}), 295) == (0.5, 295, 10)
        print('selftest ok')
    else: main(*sys.argv[1:])
