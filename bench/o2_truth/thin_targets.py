#!/usr/bin/env python3
"""Thin the reads of chosen family copies and leave their relatives untouched (register 1980, part B).

    thin_targets.py FAM.copies.tsv BAM OUT.bam OUT.targets.tsv [--min-reads 20] [--keep 3] [--seed 1980] [--threads N]

Candidates = copies of FAM.copies.tsv with >= --min-reads primary alignments of mapping quality >= 1 on >= 50 of their
exon bases. A random half become targets. Reads whose primary is at a target are cut to --keep per target; a removed
read loses every record. OUT.targets.tsv: locus, chrom, start, end, exons, reads before, reads kept."""
import sys, random, argparse, collections

def main():
    import pysam
    ap = argparse.ArgumentParser(); ap.add_argument('copies'); ap.add_argument('bam'); ap.add_argument('out'); ap.add_argument('targets')
    ap.add_argument('--min-reads', type=int, default=20); ap.add_argument('--keep', type=int, default=3); ap.add_argument('--seed', type=int, default=1980); ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args(); B = pysam.AlignmentFile(a.bam, threads=a.threads); refs = set(B.references); cand = []
    for l in list(open(a.copies))[1:]:
        f = l.rstrip('\n').split('\t')
        if f[3] not in refs: continue
        ex = [tuple(map(int, x.split('-'))) for x in f[9].split(',')]; names = set()
        for r in B.fetch(f[3], ex[0][0], ex[-1][1]):
            if r.is_secondary or r.is_supplementary or r.mapping_quality < 1: continue
            if sum(max(0, min(b, y) - max(x_, x)) for x_, b in r.get_blocks() for x, y in ex) >= 50: names.add(r.query_name)
        if len(names) >= a.min_reads: cand.append((f[12], f[3], ex, sorted(names)))
    rng = random.Random(a.seed); rng.shuffle(cand); targets = cand[:len(cand) // 2]; drop = set(); kept = {}
    for g, c, ex, names in targets: kept[g] = set(rng.sample(names, a.keep))
    allkeep = set().union(*kept.values()) if kept else set()
    for g, c, ex, names in targets: drop |= set(names) - allkeep
    with open(a.targets, 'w') as w:
        w.write('locus\tchrom\tstart\tend\texons\treads_before\treads_kept\n')
        for g, c, ex, names in targets: w.write(f"{g}\t{c}\t{ex[0][0]}\t{ex[-1][1]}\t{','.join(f'{x}-{y}' for x, y in ex)}\t{len(names)}\t{len(set(names) - drop)}\n")
    n = k = 0
    with pysam.AlignmentFile(a.out, 'wb', template=B, threads=a.threads) as w:
        for r in pysam.AlignmentFile(a.bam, threads=a.threads):
            n += 1
            if r.query_name not in drop: w.write(r); k += 1
    pysam.index(a.out)
    print(f'thin_targets: {len(cand)} candidate copies, {len(targets)} targets, {len(drop)} reads removed, {k} of {n} records kept -> {a.out}')

if __name__ == '__main__': main()
