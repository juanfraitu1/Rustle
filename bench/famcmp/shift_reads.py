#!/usr/bin/env python3
"""PREREG row 1728: read-level features of the moved intron K of S1 chains against the partner's intron J.

  shift_reads.py CLASS_TSV BAM FASTA OUT_TSV [THREADS]   |   --selftest"""
import sys, csv, ast, collections, multiprocessing
CAP = 300; NEAR = 10

def introns(read):
    """-> [(start, end, indel within NEAR bp of either end)] of one aligned read (cigar N = intron, I / D = indel)."""
    pos = read.reference_start; out = []; ind = []
    for op, n in read.cigartuples:
        if op == 3: out.append((pos, pos + n))
        elif op == 1: ind.append(pos)
        elif op == 2: ind.extend((pos, pos + n))
        if op in (0, 2, 3, 7, 8): pos += n
    return [(a, b, any(a - NEAR <= x <= a or b <= x <= b + NEAR for x in ind)) for a, b in out]

def homopolymer(seq, at):
    w = seq[max(0, at - 9):at + 9].upper(); return any(b * 5 in w for b in 'ACGT')

def count(bam, contig, K, J):
    """-> (reads carrying K, carrying J, of those with an indel beside the intron: K, J); primary records, CAP each."""
    n = collections.Counter(); d = collections.Counter(); at = K[0] if K[0] != J[0] else K[1]
    for r in bam.fetch(contig, max(0, at - 1), at + 1):
        if r.flag & 2308: continue
        for a, b, x in introns(r):
            for nm, I in (('K', K), ('J', J)):
                if (a, b) == I and n[nm] < CAP: n[nm] += 1; d[nm] += x
        if n['K'] >= CAP and n['J'] >= CAP: break
    return n['K'], n['J'], d['K'], d['J']

def work(job):
    import pysam
    bam_p, fa_p, contig, pairs = job; bam = pysam.AlignmentFile(bam_p); fa = pysam.FastaFile(fa_p); res = {}
    for K, J in pairs:
        at = K[0] if K[0] != J[0] else K[1]
        res[K, J] = (*count(bam, contig, K, J), int(homopolymer(fa.fetch(contig, max(0, at - 9), at + 9), min(9, at))))
    return contig, res

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        class R: reference_start = 100; cigartuples = [(0, 50), (1, 2), (0, 5), (3, 1000), (0, 30), (2, 3), (0, 20), (3, 500), (0, 40)]
        assert introns(R) == [(155, 1155, True), (1208, 1708, False)] and homopolymer('ACGTAAAAACGTACGTAC', 9) and not homopolymer('ACGTACGTACGTACGTAC', 9)
        print('selftest ok'); sys.exit()
    cls, bam, fa, out = sys.argv[1:5]; thr = int(sys.argv[5]) if len(sys.argv) > 5 else 16; csv.field_size_limit(1 << 30)
    rows = [r for r in csv.DictReader(open(cls), delimiter='\t') if r['class'] == 'S1']; per = collections.defaultdict(set)
    for r in rows: r['K'] = ast.literal_eval(r['introns'])[0]; r['J'] = ast.literal_eval(r['partner_introns'])[0]; per[r['contig']].add((r['K'], r['J']))
    with multiprocessing.Pool(thr) as pool: F = dict(pool.imap_unordered(work, [(bam, fa, c, sorted(p)) for c, p in per.items()]))
    with open(out, 'w') as w:
        print('contig\tstrand\tshift\treads\tpartner_reads\tin_ref\tintron_annotated\tpartner_in_ref\tnK\tnJ\tindelK\tindelJ\thomopolymer', file=w)
        for r in rows: print(*(r[k] for k in ('contig', 'strand', 'shift', 'reads', 'partner_reads', 'in_ref', 'intron_annotated', 'partner_in_ref')), *F[r['contig']][r['K'], r['J']], sep='\t', file=w)
    T = [dict(zip(l, map(int, v.split('\t')[2:]))) for l in [open(out).readline().split()[2:]] for v in open(out).read().splitlines()[1:]]
    ok = [t for t in T if t['nK'] > 0]; pairs = sum(t['in_ref'] and t['partner_in_ref'] for t in T)
    print(f'S1 chains {len(T)}; with a K read found {len(ok)}; annotated pairs {pairs}; class in RefSeq {sum(t["in_ref"] for t in T) / len(T):.3f}, intron annotated {sum(t["intron_annotated"] for t in T) / len(T):.3f}')
    f1 = lambda t: t['indelK'] / t['nK'] >= 0.5; f3 = lambda t: t['homopolymer'] == 1
    cells = [('F1 indel_K >= 0.5', f1), ('F2 indel excess >= 0.25', lambda t: t['indelK'] / t['nK'] - (t['indelJ'] / t['nJ'] if t['nJ'] else 0) >= 0.25), ('F3 homopolymer', f3),
             ('F4 intron ratio < 0.05', lambda t: t['nK'] < 0.05 * t['nJ']), ('F1 or F3', lambda t: f1(t) or f3(t))]
    print('cell\tchains\tin_RefSeq\tintron_annotated\tpairs_in_cell\tqualifies')
    for nm, f in cells:
        for sn, g in (('all shifts', lambda t: True), ('<= 10 bp', lambda t: t['shift'] <= 10)):
            v = [t for t in ok if f(t) and g(t)]; n = len(v) or 1; a = sum(t['in_ref'] for t in v) / n; b = sum(t['intron_annotated'] for t in v) / n; p = sum(t['in_ref'] and t['partner_in_ref'] for t in v) / max(1, pairs)
            print(f'{nm}, {sn}\t{len(v)}\t{a:.3f}\t{b:.3f}\t{p:.4f}\t{"yes" if len(v) >= 500 and a <= 0.033 and b <= 0.15 and p <= 0.02 else "no"}')
