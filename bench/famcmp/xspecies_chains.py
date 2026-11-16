#!/usr/bin/env python3
"""Addendum 3 of PREREG_ptr_ppa_frozen_rules_2026-11-15: chains of one sample confirmed by exact reads of a sister species.

  xspecies_chains.py MAP.paf LIFT.sh T.bam OUT.tsv Q.gtf [Q2.gtf ...]   |   --selftest
OUT: chrom, key (intron chain, 0-based half-open, as chain_rules_feat.py), lifted (1 = every intron lifts once, one contig),
sister_reads (primary reads of T with exactly the lifted chain)."""
import sys, re, collections, subprocess, pysam

def chains(paths):
    out = set()
    for p in paths:
        ex = collections.defaultdict(list)
        for l in open(p):
            f = l.split('\t')
            if len(f) > 8 and f[2] == 'exon': ex[(f[0], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
        for (c, _), e in ex.items():
            if len(e) > 1: e.sort(); out.add((c, tuple((a[1], b[0] - 1) for a, b in zip(e, e[1:]))))
    return out

def lift_chain(c, ch, L):
    """-> (contig, sorted lifted introns) or None; L = {(contig, intron): (tcontig, start, end)} of introns lifting once."""
    t = [L.get((c, i)) for i in ch]
    if None in t or len({x[0] for x in t}) != 1: return None
    return t[0][0], tuple(sorted((a, b) for _, a, b in t))

def main(paf, lift, tbam, out, *gtfs):
    C = chains(gtfs); wd = out + '.introns.bed'
    with open(wd, 'w') as o:
        for c, i in sorted({(c, i) for c, ch in C for i in ch}): o.write(f'{c}\t{i[0]}\t{i[1]}\n')
    L = {}
    for l in subprocess.run(['bash', lift, paf, wd], capture_output=True, text=True, check=True).stdout.splitlines():
        f = l.split('\t'); c, a, b = f[3].rsplit('_', 2); L[(c, (int(a), int(b)))] = (f[0], int(f[1]), int(f[2]))
    R = collections.Counter()
    for r in pysam.AlignmentFile(tbam, threads=8):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        pos = r.reference_start; ch = []
        for op, n in r.cigartuples:
            if op == 3: ch.append((pos, pos + n))
            if op in (0, 2, 3, 7, 8): pos += n
        if ch: R[(r.reference_name, tuple(ch))] += 1
    with open(out, 'w') as o:
        o.write('chrom\tkey\tlifted\tsister_reads\n')
        for c, ch in sorted(C):
            t = lift_chain(c, ch, L); o.write(f"{c}\t{','.join(f'{a}-{b}' for a, b in ch)}\t{int(t is not None)}\t{R[t] if t else 0}\n")

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        L = {('q', (10, 20)): ('t', 110, 120), ('q', (30, 40)): ('t', 130, 140), ('q', (50, 60)): ('u', 5, 15)}
        assert lift_chain('q', ((10, 20), (30, 40)), L) == ('t', ((110, 120), (130, 140)))
        assert lift_chain('q', ((10, 20), (50, 60)), L) is None and lift_chain('q', ((10, 20), (70, 80)), L) is None
        L[('q', (30, 40))] = ('t', 90, 100); assert lift_chain('q', ((10, 20), (30, 40)), L) == ('t', ((90, 100), (110, 120)))  # inversion
        print('selftest ok')
    else: main(*sys.argv[1:])
