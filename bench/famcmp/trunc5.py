#!/usr/bin/env python3
"""Opt-in 5' truncation drop (register 1677, rule 1; driver stage `trunc5`).

  trunc5.py GTF BAM OUT_GTF [THREADS]

OUT_GTF = GTF without the multi-exon transcripts whose intron chain has at least one full-chain primary read and at
least as many primary reads that run through the whole chain and continue 5' of it. Single-exon transcripts and
chains no primary read spells are kept. Last stdout line is the summary."""
import sys, re, collections, multiprocessing, pysam

def work(job):
    bam, contig, cand = job  # cand: chain -> strand
    cnt = collections.Counter()
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for op, ln in r.cigartuples:
            if op == 3: ch.append((pos, pos + ln))
            if op in (0, 2, 3, 7, 8): pos += ln
        cnt[tuple(ch)] += 1
    return contig, truncated(cnt, cand)

def truncated(cnt, cand):
    """cnt: read chain -> primary reads; cand: model chain -> strand. Chains with >= 1 full-chain read and at least as
    many reads whose chain strictly contains the chain as a contiguous run and continues on its 5' side."""
    first = collections.defaultdict(list)
    for c in cand: first[c[0]].append(c)
    cont5 = collections.Counter()
    for R, n in cnt.items():
        for p, i in enumerate(R):
            for c in first.get(i, ()):
                if len(c) < len(R) and R[p:p + len(c)] == c and (p > 0 if cand[c] == '+' else p + len(c) < len(R)): cont5[c] += n
    return {c for c in cand if cnt.get(c, 0) >= 1 and cont5[c] >= cnt[c]}

if __name__ == '__main__':
    gtf, bam, out = sys.argv[1:4]; thr = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    ex = collections.defaultdict(list)
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
    chain = {}; per = collections.defaultdict(dict)
    for k, e in ex.items():
        if len(e) < 2: continue
        e.sort(); c = tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)); chain[k] = c; per[k[0]][c] = k[1]
    have = set(pysam.AlignmentFile(bam).references)
    with multiprocessing.Pool(thr) as pool:
        drop = dict(pool.imap_unordered(work, [(bam, ctg, cand) for ctg, cand in sorted(per.items(), key=lambda x: -len(x[1])) if ctg in have]))
    gone = {k for k, c in chain.items() if c in drop.get(k[0], ())}
    with open(out, 'w') as w:
        for l in open(gtf):
            f = l.split('\t'); m = re.search(r'transcript_id "([^"]+)"', f[8]) if len(f) > 8 else None
            if not (m and (f[0], f[6], m.group(1)) in gone): w.write(l)
    print(f'dropped {len(gone)} of {len(chain)} multi-exon transcripts ({sum(len(v) for v in drop.values())} chains) as 5\' truncations')
