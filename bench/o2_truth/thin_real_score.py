#!/usr/bin/env python3
"""Truth chains recovered at thinned targets (register 1980, part B).

    thin_real_score.py TARGETS.tsv FULL_PRIMARY.gtf LABEL=GTF ...

Truth of a target = the spliced intron chains of FULL_PRIMARY.gtf (primaries-only assembly of the full reads) with >= 50
exon bases on the target. Per label: truth chains found, chains built at the targets, share of them in the truth."""
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import thin_copy_score as S

def at_targets(gtf, T):
    out = set()
    for t, (c, e, g) in S.chains(gtf).items():
        if len(e) < 2: continue
        for tc, tex in T:
            if tc == c and tex[0][0] < e[-1][1] and tex[-1][1] > e[0][0] and sum(max(0, min(y, q) - max(x, p)) for x, y in e for p, q in tex) >= 50: out.add(S.key(c, e)); break
    return out

def main(targets, full, *args):
    T = []
    for l in list(open(targets))[1:]:
        f = l.rstrip('\n').split('\t'); T.append((f[1], [tuple(map(int, x.split('-'))) for x in f[4].split(',')]))
    truth = at_targets(full, T); print(f'targets {len(T)}; truth chains (full reads, primaries only) {len(truth)}')
    print('arm\ttruth chains found\trecall\tchains built at the targets\tin the truth\tshare')
    for a in args:
        lab, p = a.split('=', 1); b = at_targets(p, T); print(f'{lab}\t{len(b & truth)}\t{len(b & truth) / max(1, len(truth)):.3f}\t{len(b)}\t{len(b & truth)}\t{len(b & truth) / max(1, len(b)):.3f}')

if __name__ == '__main__': main(*sys.argv[1:])
