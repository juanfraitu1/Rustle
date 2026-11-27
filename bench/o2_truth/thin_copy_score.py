#!/usr/bin/env python3
"""Isoforms found at the thinly planted copies only (register 1980).

    thin_copy_score.py PLANT_DIR SCENARIO LABEL=GTF ... [--thin 5]

Thin copies = planted copies with 1..--thin reads in the scenario. Per label: planted isoforms of those copies (spliced,
with >= 1 read) found by exact intron chain, split by the isoform's own planted reads (>= 2, 1); transcripts built at
those copies (>= 50 exon bases, spliced) and the share of them that are planted isoforms of that scenario."""
import sys, re, collections

def chains(path):
    ex = collections.defaultdict(list); meta = {}
    for l in open(path):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); ex[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], re.search(r'gene_id "([^"]+)"', f[8]).group(1))
    return {t: (meta[t][0], sorted(e), meta[t][1]) for t, e in ex.items()}

def key(c, e): return (c, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)))

def main(plant, sc, *args):
    thin = 5; args = list(args)
    if '--thin' in args: i = args.index('--thin'); thin = int(args[i + 1]); del args[i:i + 2]
    plan = [l.rstrip('\n').split('\t') for l in open(plant + '/plan.tsv')]; h = plan[0].index(sc); per = collections.Counter(); treads = {}
    for r in plan[1:]: per[r[1]] += int(r[h]); treads[r[3]] = int(r[h])
    T = chains(plant + '/loci.gtf'); thin_copies = {g for g, n in per.items() if 0 < n <= thin}
    truth = {}; copy_ex = collections.defaultdict(list)
    for t, (c, e, g) in T.items():
        copy_ex[g] += e
        if g in thin_copies and len(e) > 1 and treads.get(t, 0) > 0: truth[key(c, e)] = treads[t]
    loc = {g: (T[next(t for t in T if T[t][2] == g)][0], sorted(set(copy_ex[g]))) for g in thin_copies}
    allk = {key(c, e) for t, (c, e, g) in T.items() if len(e) > 1 and treads.get(t, 0) > 0}
    n2 = sum(1 for v in truth.values() if v >= 2); n1 = sum(1 for v in truth.values() if v == 1)
    print(f'## {sc} {plant}: thin copies {len(thin_copies)}, their planted spliced isoforms {len(truth)} ({n2} with 2 or more reads, {n1} with 1 read)')
    print('arm\tisoforms found\trecall\tof those with >= 2 reads\tof those with 1 read\ttranscripts at thin copies\tplanted isoforms among them\tshare')
    for a in args:
        lab, p = a.split('=', 1); A = chains(p); at = []
        for t, (c, e, g) in A.items():
            if len(e) < 2: continue
            if any(lc == c and le[0][0] < e[-1][1] and le[-1][1] > e[0][0] and sum(max(0, min(y, q) - max(x, p_)) for x, y in e for p_, q in le) >= 50 for lc, le in loc.values()): at.append(key(c, e))
        got = set(at) & set(truth); f2 = sum(1 for k in got if truth[k] >= 2); f1 = sum(1 for k in got if truth[k] == 1); good = sum(1 for k in set(at) if k in allk)
        print(f'{lab}\t{len(got)}\t{len(got) / max(1, len(truth)):.3f}\t{f2} ({f2 / max(1, n2):.3f})\t{f1} ({f1 / max(1, n1):.3f})\t{len(set(at))}\t{good}\t{good / max(1, len(set(at))):.3f}')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: assert key('c', [(0, 10), (20, 30), (40, 50)]) == ('c', ((10, 20), (30, 40))); print('thin_copy_score selftest ok')
    else: main(*sys.argv[1:])
