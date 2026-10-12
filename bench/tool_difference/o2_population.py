#!/usr/bin/env python3
"""How many multi-placement reads O2 can assign, and how many it moves off the aligner's primary (docs/PREREG_o2_population_2026-10-12.md).

usage: o2_population.py --bam B --fasta F --molecules M.tsv --n 20000 --seed 1 --out PREFIX"""
import argparse
import collections
import math
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from s1_o2 import COMP, aligned_columns, certify_pair, parse_cigar  # noqa: E402


def molecule_class(pairs):
    """pairs = [(verdict 'A' | 'B' | None, distinguishing columns)] of the primary vs each rival -> class"""
    if any(v == 'B' for v, _ in pairs):
        return 'moved'
    if pairs and all(v == 'A' for v, _ in pairs):
        return 'confirmed'
    if all(m == 0 for _, m in pairs):
        return 'identical'
    return 'undecided'


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 1.0
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def main():
    import pysam
    ap = argparse.ArgumentParser()
    for a in ('--bam', '--fasta', '--molecules', '--out'):
        ap.add_argument(a, required=True)
    ap.add_argument('--n', type=int, default=20000); ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()
    pop, total = [], 0
    for ln in open(args.molecules):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t'); total += 1
        best, second = int(f[1]), int(f[2])
        if int(f[3]) >= 2 and second >= 0.98 * best:
            pop.append((f[0], best, second, int(f[5])))
    sample = random.Random(args.seed).sample(pop, min(args.n, len(pop)))
    info = {s[0]: s for s in sample}
    nf = args.out + '.names.txt'
    with open(nf, 'w') as o:
        o.write('\n'.join(sorted(info)) + '\n')
    recs = collections.defaultdict(list)
    p = subprocess.Popen(['samtools', 'view', '-F', '2052', '-N', nf, args.bam], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for ln in p.stdout:
        f = ln.split('\t', 11)
        i = ln.find('\tAS:i:')
        a = int(ln[i + 6:].split('\t', 1)[0]) if i >= 0 else None
        recs[f[0]].append((int(f[1]), f[2], int(f[3]) - 1, parse_cigar(f[5]), a))
    p.wait()
    fa = pysam.FastaFile(args.fasta)

    def bases(chrom, cols, reverse):
        lo = min(r for r, _ in cols.values()); hi = max(r for r, _ in cols.values()) + 1
        seq = fa.fetch(chrom, lo, hi).upper()
        return (lambda o: seq[cols[o][0] - lo].translate(COMP)) if reverse else (lambda o: seq[cols[o][0] - lo])

    out = []
    for name, (_, best, second, pas) in info.items():
        rs = recs.get(name, [])
        prim = [r for r in rs if not r[0] & 256]
        rivals = [r for r in rs if r[0] & 256 and r[4] is not None and r[4] >= 0.98 * best]
        if not prim or not rivals:
            out.append((name, 'no_rival_record', pas, best, second, 0, 0, 0)); continue
        pf, pc, pp, pcig, _ = prim[0]
        ca = aligned_columns(pp, pcig, bool(pf & 16)); ba = bases(pc, ca, bool(pf & 16)) if ca else None
        pairs, mosaic, spliced_win = [], 0, 0
        for rf, rc, rp, rcig, ra in rivals:
            cb = aligned_columns(rp, rcig, bool(rf & 16))
            if not ca or not cb:
                pairs.append((None, 0)); continue
            v, ka, kb, m, _ = certify_pair(ca, cb, ba, bases(rc, cb, bool(rf & 16)))
            pairs.append((v, m)); mosaic |= ka >= 2 and kb >= 2
            if v == 'B':
                spliced_win |= sum(1 for op, _ in rcig if op == 3) > sum(1 for op, _ in pcig if op == 3)
        out.append((name, molecule_class(pairs), pas, best, second, len(rivals), int(mosaic), int(spliced_win)))
    with open(args.out + '.molecules.tsv', 'w') as o:
        o.write('read\tclass\tprimary_as\tbest_as\tsecond_as\trivals\tmosaic\tmoved_to_more_spliced\n')
        for r in out:
            o.write('\t'.join(map(str, r)) + '\n')
    scored = [r for r in out if r[1] != 'no_rival_record']; n = len(scored)
    print(f"molecules {total:,}; population (second-best AS >= 0.98 x best, >= 2 records) {len(pop):,} ({len(pop) / total:.1%}); "
          f"sampled {len(sample):,}, scored {n:,} (no rival record {len(out) - n:,})")
    c = collections.Counter(r[1] for r in scored)
    for k in ('moved', 'confirmed', 'identical', 'undecided'):
        lo, hi = wilson(c[k], n)
        print(f"  {k:10s} {c[k]:6,}  {c[k] / max(1, n):6.1%}  [{lo:.1%}, {hi:.1%}]  ~{c[k] / max(1, n) * len(pop):,.0f} in the population")
    rank = lambda r: 'true tie' if r[2] == r[3] == r[4] else ('primary below best' if r[2] < r[3] else 'primary best, close rival')
    for g in ('true tie', 'primary below best', 'primary best, close rival'):
        sub = [r for r in scored if rank(r) == g]; cc = collections.Counter(r[1] for r in sub)
        print(f"  {g:26s} {len(sub):6,}: " + ', '.join(f"{k} {cc[k] / max(1, len(sub)):.1%}" for k in ('moved', 'confirmed', 'identical', 'undecided')))
    mv = [r for r in scored if r[1] == 'moved']
    print(f"  mosaic {sum(r[6] for r in scored) / max(1, n):.1%}; moved molecules whose winning placement has more introns than the primary: "
          f"{sum(r[7] for r in mv):,} of {len(mv):,}")


if __name__ == '__main__':
    main()
