#!/usr/bin/env python3
"""Tag multi-exon transcripts whose longest ORF stops early (register 1766, rule P1 as an attribute).

    gtf_ptc.py IN.gtf GENOME.fa OUT.gtf   |   gtf_ptc.py --selftest

Appends to the transcript line of every stranded multi-exon transcript:
  orf_len "N"   length in nt of the longest ATG..stop ORF on the transcript strand (0 = none)
  ptc "1"       that ORF's stop lies more than 50 nt before the last exon junction, else "0"
Nothing is added, dropped or reordered. On seven ape samples chains with ptc "1" match RefSeq 3-5 times less often than
the rest, but about half of the correct non-coding (NR / XR) chains carry it: a tag to filter on, not a verdict.
Prints one summary line (the driver logs it)."""
import collections, re, sys

COMP = str.maketrans('ACGTN', 'TGCAN'); STOP = {'TAA', 'TAG', 'TGA'}

def orf(s):
    """-> (length incl. stop, end offset) of the longest ATG..stop in s"""
    best = (0, 0)
    for fr in range(3):
        start = None
        for i in range(fr, len(s) - 2, 3):
            c = s[i:i + 3]
            if start is None and c == 'ATG': start = i
            elif start is not None and c in STOP:
                if i + 3 - start > best[0]: best = (i + 3 - start, i + 3)
                start = None
    return best

def call(exons, strand, fetch):
    """exons 0-based half-open, sorted; fetch(start, end) -> sequence. -> (orf_len, ptc)"""
    s = ''.join(fetch(a, b) for a, b in exons).upper()
    last_exon = exons[-1] if strand == '+' else exons[0]
    if strand == '-': s = s.translate(COMP)[::-1]
    n, end = orf(s)
    return n, int(n > 0 and end < len(s) - (last_exon[1] - last_exon[0]) - 50)

def main(gtf, fasta, out):
    import pysam
    fa = pysam.FastaFile(fasta); lines = open(gtf).read().splitlines(True); ex = collections.defaultdict(list)
    tid = lambda a: re.search(r'transcript_id "([^"]+)"', a).group(1)
    for l in lines:
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[tid(f[8])].append((int(f[3]) - 1, int(f[4])))
    n = k = 0
    with open(out, 'w') as w:
        for l in lines:
            f = l.rstrip('\n').split('\t')
            if len(f) > 8 and f[2] == 'transcript' and f[6] in '+-' and len(ex[tid(f[8])]) > 1:
                o, p = call(sorted(ex[tid(f[8])]), f[6], lambda a, b, c=f[0]: fa.fetch(c, a, b)); n += 1; k += p
                l = l.rstrip('\n').rstrip() + ('' if l.rstrip().endswith(';') else ';') + ' orf_len "%d"; ptc "%d";\n' % (o, p)
            w.write(l)
    print('%d of %d multi-exon transcripts tagged ptc "1"' % (k, n))

def selftest():
    assert orf('CCATGAAATAGCC') == (9, 11) and orf('CCC') == (0, 0)
    # exon 1 = ATG AAA TAG + 60 filler, exon 2 = 10 nt: the stop is 60 nt before the junction
    g = 'ATGAAATAG' + 'C' * 60 + 'NNNNN' + 'C' * 10; f = lambda a, b: g[a:b]; e = [(0, 69), (74, 84)]
    assert call(e, '+', f) == (9, 1)
    # the same transcript on the minus strand of the reverse complement
    r = g.translate(COMP)[::-1]; fr = lambda a, b: r[a:b]; er = [(0, 10), (15, 84)]
    assert call(er, '-', fr) == (9, 1)
    # stop in the last exon: no flag
    g2 = 'ATGAAA' + 'NNNNN' + 'AAATAGCC'; assert call([(0, 6), (11, 19)], '+', lambda a, b: g2[a:b]) == (12, 0)
    print('gtf_ptc selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    elif len(sys.argv) == 4: main(*sys.argv[1:])
    else: sys.exit(__doc__)
