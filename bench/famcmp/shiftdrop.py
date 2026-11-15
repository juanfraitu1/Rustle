#!/usr/bin/env python3
"""Opt-in drop of weaker shifted-site chains with indels beside the moved intron (registers 1728-1730, 1732, 1734, 1738, 1739;
driver stage `shiftdrop`).

  shiftdrop.py GTF BAM OUT_GTF [THREADS]   |   --selftest

A chain is "weaker shifted-site" when another chain of the GTF on the same contig and strand has at least its reads
and is equal to it except one intron that shares exactly one boundary with its own (one internal exon end moved);
terminal-exon ends are never touched. OUT_GTF = GTF without the transcripts of such chains for which at least half of
the primary reads carrying the moved intron have an insertion or deletion within 10 bp of it.
Such chains are in RefSeq 2-4% of the time and their moved intron is annotated (RefSeq + CAT) 7-10% of the time on two
gorilla and two human samples. It is a small clean-up of the isoform catalogue (0.3-0.9% of chains); no other
stage reads its output. Uses no annotation. Last stdout line is the summary."""
import sys, re, collections, multiprocessing
from alt_events import load
from shift_class import index, partners
from shift_reads import count

def indel_cell(n_reads, n_indel):
    return n_reads > 0 and n_indel / n_reads >= 0.5

def work(job):
    import pysam
    bam_p, key, d = job; bam = pysam.AlignmentFile(bam_p); one, two = index(d); gone = {}
    for c, n in d.items():
        ps = [p for p in partners(c, one, two) if d[p[2]] > n or (d[p[2]] == n and p[2] < c)]
        if not ps: continue
        cls, new, par, sh = min(ps, key=lambda p: (p[0] == 'S1t', -d[p[2]]))  # as shift_class.py
        if cls != 'S1': continue
        nK, _, dK, _ = count(bam, key[0], new[0], sorted(set(par) - set(c))[0])
        if indel_cell(nK, dK): gone[c] = sh
    return key, gone

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        assert indel_cell(2, 1) and indel_cell(14, 14) and not indel_cell(0, 0) and not indel_cell(7, 3)
        print('selftest ok'); sys.exit()
    gtf, bam, out = sys.argv[1:4]; thr = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    D = load(gtf, 'reads'); jobs = [(bam, k, {c: max(t[2] or 0 for t in v) for c, v in D[k].items()}) for k in sorted(D)]
    with multiprocessing.Pool(thr) as pool: G = dict(pool.imap_unordered(work, jobs))
    ex = collections.defaultdict(list); T = re.compile(r'transcript_id "([^"]+)"'); L = open(gtf).readlines()
    for l in L:
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[f[0], f[6], T.search(f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
    bad = {k[2] for k, e in ex.items() if (e.sort() or True) and tuple((a[1], b[0]) for a, b in zip(e, e[1:])) in G.get(k[:2], ())}
    with open(out, 'w') as w: w.writelines(l for l in L if not ((m := T.search(l)) and m.group(1) in bad))
    n = sum(len(g) for g in G.values())
    print(f"{n} of {sum(len(v) for v in D.values())} multi-exon chains dropped ({len(bad)} transcripts): weaker shifted-site chains with an indel beside the moved intron in >= half of its reads")
