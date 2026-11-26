#!/usr/bin/env python3
"""How much an assembly changes between two runs (register 1969).

    stability.py PLANT_DIR SCENARIO LABEL=A.gtf,B.gtf ... [--expressed LABEL=A.o2.gtf,B.o2.gtf]

Per label: spliced intron chains only in A / only in B / in both; planted copies with a transcript (>= 50 exon bases of
the copy) only in A / only in B / in both; the same for the SILENT copies of the scenario (0 planted reads).
--expressed: the copies holding a transcript with o2_status "expressed" instead of every transcript."""
import sys, os, re, collections

def load(path, only=None):
    ex = collections.defaultdict(list); meta = {}; ok = set()
    for l in open(path):
        f = l.split('\t')
        if len(f) < 9: continue
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if not t: continue
        t = t.group(1)
        if f[2] == 'transcript' and (only is None or only in f[8]): ok.add(t)
        if f[2] == 'exon': ex[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6])
    return {t: (meta[t][0], meta[t][1], sorted(e)) for t, e in ex.items() if only is None or t in ok}

def chains(tx): return {(c, s, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))) for c, s, e in tx.values() if len(e) > 1}

def copies(tx, loci):
    out = set()
    for c, s, e in tx.values():
        for g, (gc, ge) in loci.items():
            if gc == c and ge[0][0] < e[-1][1] and ge[-1][1] > e[0][0] and sum(max(0, min(b, y) - max(a, x)) for a, b in e for x, y in ge) >= 50: out.add(g)
    return out

def main(plant, sc, *args):
    plan = [l.rstrip('\n').split('\t') for l in open(plant + '/plan.tsv')]; h = plan[0].index(sc); n = collections.Counter()
    for r in plan[1:]: n[r[1]] += int(r[h])
    silent = {g for g, v in n.items() if v == 0}; ex = collections.defaultdict(list); lc = {}
    for l in open(plant + '/loci.gtf'):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); ex[g].append((int(f[3]) - 1, int(f[4]))); lc[g] = f[0]
    loci = {g: (lc[g], sorted(set(e))) for g, e in ex.items()}
    print('tool\tchains: only first / only second / both\tcopies with a transcript: only first / only second / both\tsilent copies called: only first / only second / both')
    only = None
    for a in args:
        if a == '--expressed': only = 'o2_status "expressed"'; continue
        lab, p = a.split('=', 1); pa, pb = p.split(',')
        if not (os.path.exists(pa) and os.path.exists(pb)): print(f'{lab}\tMISSING'); continue
        A, B = load(pa, only), load(pb, only); ca, cb = chains(A), chains(B); ga, gb = copies(A, loci), copies(B, loci)
        print(f"{lab}{' (labelled expressed)' if only else ''}\t{len(ca - cb)} / {len(cb - ca)} / {len(ca & cb)}\t{len(ga - gb)} / {len(gb - ga)} / {len(ga & gb)}\t{len((ga - gb) & silent)} / {len((gb - ga) & silent)} / {len(ga & gb & silent)}")
        only = None

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        t = {'a': ('c', '+', [(0, 100), (200, 300)]), 'b': ('c', '+', [(1000, 1100)])}
        assert chains(t) == {('c', '+', ((100, 200),))} and copies(t, {'G': ('c', [(50, 250)]), 'H': ('c', [(1090, 1200)]), 'I': ('d', [(0, 100)])}) == {'G'}; print('stability selftest ok')
    else: main(*sys.argv[1:])
