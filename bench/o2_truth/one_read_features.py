#!/usr/bin/env python3
"""Features of the lone read behind each one-read locus, with its labels (register 1987 addendum).

    one_read_features.py ONELOCI.gtf BASE.gtf BAM GENOME.fa TRUTH.gtf REP1.bam[,REP2.bam] OUT.tsv [THREADS]

ONELOCI.gtf = PREFIX.oneloci.gtf (transcripts tagged one_read_locus). One line per tagged transcript: the lone primary read's
features and whether the chain is annotated / seen in a replicate."""
import sys, os, bisect, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread, one_read_loci as O
RC = str.maketrans('ACGT', 'TGCA')

def main(gtf, base, bam, fa, truth, reps, out, thr=32):
    import pysam
    M = oneread.models(gtf); C = {t: v for t, v in M.items() if any('one_read_locus' in l for l in v[3])}
    keys = {(c, oneread.chain(ex)) for c, s, ex, _ in C.values()}
    tr = {(c, oneread.chain(ex)) for c, s, ex, _ in oneread.models(truth).values() if len(ex) > 1}
    rep = set()
    for b in reps.split(','): rep |= set(O.reads_with(b, keys, int(thr)))
    sp = collections.defaultdict(list)
    for t, (c, s, ex, _) in M.items():
        if t not in C: sp[c].append((ex[0][0], ex[-1][1], s))
    for v in sp.values(): v.sort()
    B = pysam.AlignmentFile(bam); F = pysam.FastaFile(fa)
    with open(out, 'w') as o:
        o.write('tid\tannotated\treplicated\tintrons\tmin_intron\tmax_intron\tmin_end_exon\tclip\tde\tmapq\tcanonical\tdist\tantisense\taligned\n')
        for t, (c, s, ex, _) in C.items():
            ch = oneread.chain(ex); r = next((a for a in B.fetch(c, ex[0][1] - 1, ex[0][1]) if not a.is_secondary and not a.is_supplementary and O.read_chain(a) == ch), None)
            if r is None: continue
            mot = [(F.fetch(c, a, a + 2) + F.fetch(c, b - 2, b)).upper() for a, b in ch]
            mot = [m if s == '+' else m.translate(RC)[::-1] for m in mot]
            ct = r.cigartuples; clip = sum(n for op, n in (ct[0], ct[-1]) if op == 4)
            L = sp[c]; i = bisect.bisect_left(L, (ex[0][0],)); near = L[max(0, i - 3):i + 3]; a0, b0 = ex[0][0], ex[-1][1]
            dist = min((max(0, x - b0, a0 - y) for x, y, _ in near), default=10**9); anti = int(any(x < b0 and a0 < y and z != s for x, y, z in near))
            o.write('\t'.join(map(str, (t, int((c, ch) in tr), int((c, ch) in rep), len(ch), min(b - a for a, b in ch), max(b - a for a, b in ch),
                    min(ch[0][0] - r.reference_start, r.reference_end - ch[-1][1]), clip, r.get_tag('de') if r.has_tag('de') else 'NA', r.mapping_quality,
                    int(all(m in ('GTAG', 'GCAG', 'ATAC') for m in mot)), dist, anti, r.query_alignment_length))) + '\n')

if __name__ == '__main__': main(*sys.argv[1:])
