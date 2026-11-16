"""orf.py GTF FASTA OUT : per multi-exon transcript, longest ATG..stop ORF on its strand and whether the stop lies > 50 nt before the last junction."""
import sys, re, collections, pysam
fa = pysam.FastaFile(sys.argv[2]); ex = collections.defaultdict(list); rd = {}
for l in open(sys.argv[1]):
    f = l.split('\t')
    if len(f) < 9: continue
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    if f[2] == 'exon': ex[(f[0], f[6], t)].append((int(f[3]) - 1, int(f[4])))
    elif f[2] == 'transcript': m = re.search(r'reads "([^"]+)"', f[8]); rd[t] = m.group(1) if m else '0'
comp = str.maketrans('ACGTN', 'TGCAN'); STOP = {'TAA', 'TAG', 'TGA'}
def orf(s):  # (length incl. stop, end offset) of the longest ATG..stop
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
with open(sys.argv[3], 'w') as w:
    w.write('chrom\tchain\ttid\treads\ttx_len\torf_len\tptc\n')
    for (c, st, t), e in ex.items():
        e.sort()
        if len(e) < 2: continue
        s = ''.join(fa.fetch(c, a, b) for a, b in e).upper()
        last = len(s) - (e[-1][1] - e[-1][0])  # offset of the last junction on +
        if st == '-': s = s.translate(comp)[::-1]; last = len(s) - (e[0][1] - e[0][0])
        n, end = orf(s)
        w.write('\t'.join(map(str, [c, ','.join('%d-%d' % (e[k][1], e[k + 1][0]) for k in range(len(e) - 1)), t, rd.get(t, 0), len(s), n, int(n > 0 and end < last - 50)])) + '\n')
