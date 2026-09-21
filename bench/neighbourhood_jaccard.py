#!/usr/bin/env python3
"""Neighbourhood Jaccard as a connectivity metric, WITH the component-size control.
Per `docs/PREREG_neighbourhood_jaccard_2026-09-21.md`.

    J_N(u,v) = |N(u) ∩ N(v)| / |N(u) ∪ N(v)|     over graph neighbours, u and v themselves excluded

r459 measured this at standalone AUC 0.826 — the highest separator in the register — but never
size-residualised it, and r523 showed edge betweenness collapses 0.683 -> 0.531 under exactly that
control (corr -0.71 with log component size). So the raw AUC is reported for comparability and the
STRATIFIED AUC is the number that decides.

r459's structural warning is also quantified here: a 2-copy family is an isolated edge with no common
neighbours, so J_N is identically zero on the modal family.

Usage: neighbourhood_jaccard.py --graphs DIR --gff G --soto S1C --chroms chr2,chr8,chr10
"""
import argparse
import collections
import csv
import itertools
import math
import re


def gene_names(gff, chrom):
    n = {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            n[f'{f[0]}:{f[3]}-{f[4]}'] = m.group(1)
    return n


def soto(s1c):
    ids = collections.defaultdict(set)
    for r in csv.DictReader(open(s1c), delimiter='\t'):
        f = (r.get('Family ID') or '').strip(); n = (r.get('Gene Name') or '').strip()
        if f and f != 'N/A' and n:
            ids[n].add(f)
    return {n: next(iter(v)) for n, v in ids.items() if len(v) == 1}


def auc(tv, fv):
    if not tv or not fv:
        return None
    pr = list(itertools.product(tv, fv))
    return (sum(1 for x, y in pr if x > y) + 0.5 * sum(1 for x, y in pr if x == y)) / len(pr)


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--soto', '--chroms'):
        ap.add_argument(x, required=True)
    ap.add_argument('--graphs', help='dir of <chrom>.graph.tsv (the SHIPPED post-conjunct graph)')
    ap.add_argument('--pafs', help='chrom=path,... : build the PRE-CONJUNCT graph from raw alignments '
                                   '(identity >= 0.7, cov_longer >= 0.3, >= 300 bp) -- the population a '
                                   'connectivity metric would actually be used on')
    a = ap.parse_args()
    fam = soto(a.soto)

    rows = []            # (label, J_N, component_size)
    famsize = collections.Counter()
    for chrom in a.chroms.split(','):
        names = gene_names(a.gff, chrom)
        adj = collections.defaultdict(set)
        edges = []
        if a.pafs:
            paf = dict(x.split('=') for x in a.pafs.split(','))[chrom]
            best = {}
            for line in open(paf):
                f = line.rstrip('\n').split('\t')
                if len(f) < 11 or f[0] == f[5]:
                    continue
                nm, al = int(f[9]), int(f[10])
                if nm < 300 or not al or nm / al < 0.7:
                    continue
                if nm / max(int(f[1]), int(f[6])) < 0.3:
                    continue
                k = tuple(sorted((f[0], f[5])))
                best[k] = max(best.get(k, 0), nm)
            for (u, v) in best:
                adj[u].add(v); adj[v].add(u); edges.append((u, v))
        else:
            for line in open(f'{a.graphs}/{chrom}.graph.tsv'):
                f = line.rstrip('\n').split('\t')
                if len(f) != 3 or f[0] == f[1]:
                    continue
                adj[f[0]].add(f[1]); adj[f[1]].add(f[0])
                edges.append((f[0], f[1]))
        # component sizes
        seen, comp = set(), {}
        for n0 in adj:
            if n0 in seen:
                continue
            stack, members = [n0], []
            while stack:
                x = stack.pop()
                if x in seen:
                    continue
                seen.add(x); members.append(x); stack.extend(adj[x] - seen)
            for m in members:
                comp[m] = len(members)
        # truth family sizes on this chromosome (for the pair fraction)
        onchrom = collections.Counter()
        for sp, g in names.items():
            if g in fam:
                onchrom[fam[g]] += 1
        for f_, n in onchrom.items():
            if n >= 2:
                famsize[n] += 1
        for u, v in edges:
            gu, gv = names.get(u), names.get(v)
            if not gu or not gv:
                continue
            fu, fv_ = fam.get(gu), fam.get(gv)
            if not fu or not fv_:
                continue
            nu, nv = adj[u] - {v}, adj[v] - {u}
            uni = nu | nv
            j = len(nu & nv) / len(uni) if uni else 0.0
            rows.append((fu == fv_, j, comp.get(u, 1)))

    T = [r for r in rows if r[0]]; F = [r for r in rows if not r[0]]
    npair = famsize[2]; nfam = sum(famsize.values())
    print(f"truth families with >=2 members: {nfam} | exactly 2 members: {npair} "
          f"({100*npair/nfam if nfam else 0:.1f}%)   [r459 cited 57%]\n")
    print(f"scored edges: {len(rows)} — {len(T)} TRUE, {len(F)} FALSE")
    cov = sum(1 for r in T if r[1] > 0)
    print(f"TRUE edges with J_N > 0 (metric not blind): {cov}/{len(T)} = {100*cov/len(T) if T else 0:.1f}%\n")
    raw = auc([r[1] for r in T], [r[1] for r in F])
    print(f"  RAW AUC                     {raw:.3f}    [r459 reported 0.826]")

    bins = [(2, 2), (3, 4), (5, 9), (10, 10**9)]
    tot_w = 0; acc = 0
    print(f"\n  {'component size':>15} {'TRUE':>6} {'FALSE':>6} {'AUC':>8}")
    for lo, hi in bins:
        tv = [r[1] for r in T if lo <= r[2] <= hi]
        fv = [r[1] for r in F if lo <= r[2] <= hi]
        A = auc(tv, fv)
        lbl = f'{lo}' if lo == hi else (f'>={lo}' if hi > 10**8 else f'{lo}-{hi}')
        print(f"  {lbl:>15} {len(tv):>6} {len(fv):>6} {'  n/a' if A is None else format(A, '>8.3f')}")
        if A is not None:
            w = len(tv) + len(fv); acc += A * w; tot_w += w
    print(f"\n  SIZE-RESIDUALISED AUC       {acc/tot_w if tot_w else 0:.3f}    <- the deciding number")
    xs = [math.log(max(r[2], 1)) for r in rows]; ys = [r[1] for r in rows]
    mx, my = sum(xs)/len(xs), sum(ys)/len(ys)
    num = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    den = math.sqrt(sum((x-mx)**2 for x in xs) * sum((y-my)**2 for y in ys))
    print(f"  corr(J_N, log component size) {num/den if den else 0:+.3f}    [r523's betweenness was -0.71]")


if __name__ == '__main__':
    main()
