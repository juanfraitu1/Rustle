#!/usr/bin/env python3
"""Splice-site strength from the sample's own introns (the PWM code of register row 1789, `ssfeat.py`, moved here).

    sites.py --selftest
Donor = 9-mer (3 exonic + 6 intronic bases), acceptor = 23-mer (20 intronic + 3 exonic); a site's strength = the log-odds
of its bases under the position weights of every intron the sample's assembly emits, against a flat background."""
import collections, math, re, sys
COMP = str.maketrans('ACGTNacgtn', 'TGCANTGCAN'); rc = lambda s: s.translate(COMP)[::-1]
PAD = 40   # bases fetched on each side of an intron end; an end closer than this to a contig end has no window

def windows(seq, a, b, strand):
    """transcript-oriented windows around intron [a, b): donor window (intron's first base at index 13), acceptor window (first exonic base at index 30)"""
    if strand == '+': return seq[a - 13:a + 16].upper(), seq[b - 30:b + 13].upper()
    return rc(seq[b - 16:b + 13]).upper(), rc(seq[a - 13:a + 30]).upper()
don = lambda w, k=13: w[k - 3:k + 6]; acc = lambda w, k=30: w[k - 20:k + 3]
def pwm(sites, L):
    c = [collections.Counter() for _ in range(L)]
    for s in sites:
        if len(s) == L:
            for i, x in enumerate(s): c[i][x] += 1
    n = sum(c[0].values()) or 1
    return [{x: math.log2((c[i][x] + 1) / (n + 4) / 0.25) for x in 'ACGT'} for i in range(L)]
score = lambda m, s: sum(m[i].get(x, -2.0) for i, x in enumerate(s)) if len(s) == len(m) else float('nan')

def end_windows(fa, intron):
    """windows() of an intron (chrom, strand, start, end; 1-based closed) read from an indexed FASTA; None near a contig end"""
    chrom, strand, s, e = intron; a, b = s - 1, e
    if a - PAD < 0 or b + PAD > fa.get_reference_length(chrom) or b <= a: return None
    left, right = fa.fetch(chrom, a - PAD, a + PAD), fa.fetch(chrom, b - PAD, b + PAD)
    return windows(left + right, PAD, 3 * PAD, strand)      # the two pieces side by side: the intron is [PAD, 3 PAD) of it

def other_end(fa, intron, anchor):
    """the site at the end the intron does NOT share with its bunch (anchor = the shared end): its 9-mer or 23-mer, or None"""
    w = end_windows(fa, intron)
    if w is None: return None
    s = acc(w[1]) if anchor == 'donor' else don(w[0])
    return None if 'N' in s else s
canonical = lambda site, anchor: site is not None and (site[18:20] == 'AG' if anchor == 'donor' else site[3:5] in ('GT', 'GC'))

def other_end_strength(fa, pwms, intron, anchor):
    s = other_end(fa, intron, anchor)
    return float('nan') if s is None else score(pwms[1] if anchor == 'donor' else pwms[0], s)

def sample_pwms(gtf, fasta):
    """(donor PWM, acceptor PWM) from the introns of the GTF's stranded multi-exon transcripts"""
    import pysam
    ex = collections.defaultdict(list); st = {}
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            t = (f[0], re.search(r'transcript_id "([^"]+)"', f[8]).group(1)); ex[t].append((int(f[3]), int(f[4]))); st[t] = f[6]
    introns = set()
    for t, e in ex.items():
        e.sort()
        if st[t] in '+-': introns.update((t[0], st[t], e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1))
    fa = pysam.FastaFile(fasta); W = [w for w in (end_windows(fa, i) for i in sorted(introns)) if w]
    return pwm([don(w[0]) for w in W], 9), pwm([acc(w[1]) for w in W], 23)

def selftest():
    g = 'A' * 40 + 'CAG' + 'GTAAGT' + 'C' * 30 + 'TTTTTTTTTTCTCTTTTCAG' + 'GAT' + 'A' * 40   # intron [43, 99)
    wd, wa = windows(g, 43, 99, '+'); assert don(wd) == 'CAGGTAAGT', don(wd); assert acc(wa) == 'TTTTTTTTTTCTCTTTTCAGGAT', acc(wa)
    r = rc(g); n = len(g); wd2, wa2 = windows(r, n - 99, n - 43, '-'); assert (wd2, wa2) == (wd, wa)
    md = pwm(['CAGGTAAGT'] * 50 + ['AAGGTGAGT'] * 50, 9); ma = pwm([acc(wa)] * 100, 23)
    assert score(md, 'CAGGTAAGT') > score(md, 'CAGATAAGT') and math.isnan(score(md, 'CAG'))
    class F:      # an indexed FASTA of one contig
        def get_reference_length(self, c): return len(g)
        def fetch(self, c, a, b): return g[a:b]
    assert end_windows(F(), ('c', '+', 44, 99)) == (wd, wa) and end_windows(F(), ('c', '+', 30, 99)) is None
    assert other_end(F(), ('c', '+', 44, 99), 'donor') == acc(wa) and other_end(F(), ('c', '+', 44, 99), 'acceptor') == don(wd)
    assert canonical(acc(wa), 'donor') and canonical(don(wd), 'acceptor') and not canonical('TTTTTTTTTTCTCTTTTCAAGAT', 'donor')
    print('selftest ok')
if __name__ == '__main__': selftest()
