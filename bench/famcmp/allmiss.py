#!/usr/bin/env python3
"""One row per RefSeq multi-exon intron chain: which arms output it and what the reads show (register 1766).

  allmiss.py REF_GTF BAM OUT_PREFIX NAME=GTF [NAME=GTF ...]
Writes OUT_PREFIX.tsv and OUT_PREFIX.reads.pkl (Counter of (chrom, chain) over primary records, reused by the rule scripts).
Read evidence ignores strand (a read carries none). Tiers:
  full     >= 1 primary record whose whole intron chain equals the chain
  linked   no such record, every intron in a primary record and every adjacent intron pair together in one
  unlinked every intron seen, some adjacent pair never in one record
  nodata   some intron in no primary record
"""
import sys, re, gzip, pickle, collections, pysam

def chains(path):
    ex = collections.defaultdict(list)
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    out = set()
    for (c, s, _), e in ex.items():
        e.sort()
        if len(e) > 1: out.add((c, s, tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))))  # intron = 0-based half-open
    return out

def read_chain(a):
    p = a.reference_start; ch = []
    for op, n in a.cigartuples:
        if op == 3: ch.append((p, p + n))
        if op in (0, 2, 3, 7, 8): p += n
    return tuple(ch)

def scan(bam):
    prim = collections.Counter(); sec = collections.Counter()
    for a in pysam.AlignmentFile(bam):
        if a.is_unmapped or a.is_supplementary: continue
        ch = read_chain(a)
        if ch: (sec if a.is_secondary else prim)[(a.reference_name, ch)] += 1
    return prim, sec

if __name__ == '__main__':
    ref_gtf, bam, pre = sys.argv[1:4]; arms = [a.split('=', 1) for a in sys.argv[4:]]
    ref = chains(ref_gtf); got = [(n, chains(p)) for n, p in arms]
    prim, sec = scan(bam); pickle.dump(prim, open(pre + '.reads.pkl', 'wb'))
    intron = collections.Counter(); pair = collections.Counter(); run = collections.Counter()  # run: primary records that are a contiguous part of a ref chain
    by_first = collections.defaultdict(set)
    for c, s, ch in ref:
        for k in range(len(ch)): by_first[(c, ch[k])].add(ch)
    for (c, ch), n in prim.items():
        for i in ch: intron[(c, i)] += n
        for k in range(len(ch) - 1): pair[(c, ch[k], ch[k + 1])] += n
        for r in by_first.get((c, ch[0]), ()):
            k = r.index(ch[0])
            if r[k:k + len(ch)] == ch and len(ch) < len(r): run[(c, r)] += n
    with open(pre + '.tsv', 'w') as w:
        w.write('chrom\tstrand\tn_introns\tfull\tfull_sec\tpartial\tmin_intron\tmin_pair\ttier\t' + '\t'.join(n for n, _ in got) + '\tchain\n')
        for c, s, ch in sorted(ref):
            f = prim[(c, ch)]; mi = min(intron[(c, i)] for i in ch)
            mp = min((pair[(c, ch[k], ch[k + 1])] for k in range(len(ch) - 1)), default=mi)
            tier = 'full' if f else 'nodata' if mi == 0 else 'linked' if mp > 0 else 'unlinked'
            w.write('\t'.join(map(str, [c, s, len(ch), f, sec[(c, ch)], run[(c, ch)], mi, mp, tier] + [int((c, s, ch) in g) for _, g in got]
                                  + [','.join('%d-%d' % i for i in ch)])) + '\n')
