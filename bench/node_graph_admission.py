#!/usr/bin/env python3
"""Where do assembled de novo loci lose their place in the family graph?
Per `docs/NODE_GRAPH_ADMISSION_2026-09-21.md` (§6w2).

Every node-construction lever since §6u7 has attacked where a locus's BOUNDARIES fall. Measured end to
end on chr16 that is not where the loss is: of 1,892 loci with real homology only 864 become graph
nodes, and 99.8% of the rest are rejected by `cov_longer` alone -- identity and the 300bp floor reject
essentially nothing.

`cov_longer`'s denominator is the LONGER locus's whole span (`annotation_families.rs:407`), so an
over-merged readthrough locus inflates that denominator for every partner it aligns to and evicts
correctly-assembled loci from the graph. That collateral loss is ~2.4x larger than the over-merge itself.

⚠ `cov_longer` MUST be computed the deferred way -- union of MERGED intervals on the longer locus over
that locus's length -- because the shipped `--min-exonic-bp 1` config takes the deferred path, which
accumulates every PAF record for a pair. A best-single-record approximation overstates the effect
(78.0% vs 66.4% longer-partner; ratio 9.79x vs 3.92x). The reimplementation here lands within 2 loci of
the shipped graph (866 vs 864), which is what licenses the decomposition.

Usage:
  node_graph_admission.py --paf dn16.paf --graph dn16.graph.tsv [--gff chr16.genes.gff]
"""
import argparse
import collections
import re
import statistics

MIN_IDENT, MIN_COV, MIN_ALEN = 0.70, 0.30, 300


def merged_len(intervals):
    intervals = sorted(intervals)
    total = 0
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s > ce:
            total += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return total + ce - cs


def load(paf):
    """pair -> (intervals on the LONGER locus, summed nmatch, summed blocklen); plus locus lengths."""
    qlen = {}
    acc = collections.defaultdict(lambda: [[], 0, 0])
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 11:
            continue
        q, t = f[0], f[5]
        la, lb = int(f[1]), int(f[6])
        qlen[q] = la
        qlen[t] = lb
        if q == t:
            continue
        key = (q, t) if q <= t else (t, q)
        iv = (int(f[2]), int(f[3])) if la >= lb else (int(f[7]), int(f[8]))
        e = acc[key]
        e[0].append(iv)
        e[1] += int(f[9])
        e[2] += int(f[10])
    return qlen, acc


def gene_counter(gff):
    spans = []
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('gene', 'pseudogene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            spans.append((int(f[3]), int(f[4]), m.group(1)))
    spans.sort()

    def n_inside(name):
        m = re.match(r'^\S+:(\d+)-(\d+)$', name)
        if not m:
            return 0
        s, e = int(m.group(1)), int(m.group(2))
        c = 0
        for gs, ge, _ in spans:
            if ge < s:
                continue
            if gs > e:
                break
            ov = min(e, ge) - max(s, gs)
            if ov > 0 and ov >= 0.5 * (ge - gs + 1):
                c += 1
        return c
    return n_inside


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--paf', required=True)
    ap.add_argument('--graph', required=True)
    ap.add_argument('--gff')
    a = ap.parse_args()

    graph = set()
    for line in open(a.graph):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2 and f[0] != f[1]:
            graph.add(f[0])
            graph.add(f[1])

    qlen, acc = load(a.paf)
    passing, best = set(), {}
    for (x, y), (iv, nmatch, blocklen) in acc.items():
        if not iv or blocklen == 0:
            continue
        cov = min(merged_len(iv) / max(max(qlen[x], qlen[y]), 1), 1.0)
        ident = nmatch / blocklen
        if cov >= MIN_COV and ident >= MIN_IDENT and blocklen >= MIN_ALEN:
            passing.add(x)
            passing.add(y)
        for u, v in ((x, y), (y, x)):
            if u not in best or cov > best[u][0]:
                best[u] = (cov, ident, blocklen, v)

    homol = {x for k in acc for x in k}
    dropped = homol - passing
    print(f"assembled loci in PAF              : {len(qlen)}")
    print(f"  with any non-self homology       : {len(homol)}")
    print(f"  clearing the 3 numeric gates     : {len(passing & homol)}   (reimplementation)")
    print(f"  actually in the shipped graph    : {len(graph & homol)}   <- agreement check")
    print(f"  dropped                          : {len(dropped)} "
          f"({100*len(dropped)/len(homol):.1f}% of homologous)\n")

    reasons = collections.Counter()
    for n in dropped:
        cov, ident, blocklen, _ = best[n]
        r = []
        if cov < MIN_COV:
            r.append('cov_longer<0.30')
        if ident < MIN_IDENT:
            r.append('identity<0.70')
        if blocklen < MIN_ALEN:
            r.append('alen<300')
        reasons[' + '.join(r) if r else 'pair-level combination'] += 1
    print("rejection reason (from each locus's own best-coverage pair):")
    for k, v in reasons.most_common():
        print(f"  {v:>5} {100*v/len(dropped):>5.1f}%  {k}")

    ratios = [qlen[best[n][3]] / qlen[n] for n in dropped]
    longer = [n for n in dropped if qlen[best[n][3]] > qlen[n]]
    print(f"\ncollateral check -- is the culprit an over-merged partner?")
    print(f"  best partner LONGER than the locus : {len(longer)} "
          f"({100*len(longer)/len(dropped):.1f}%)")
    print(f"  median partner/self length ratio   : {statistics.median(ratios):.2f}")
    print(f"  median length of the evicted locus : {statistics.median([qlen[n] for n in dropped]):.0f} bp")

    if a.gff:
        n_inside = gene_counter(a.gff)
        culprit = collections.Counter(best[n][3] for n in longer)
        total = sum(culprit.values())
        multi = sum(c for name, c in culprit.items() if n_inside(name) >= 2)
        print(f"  culprit partner holds >=2 whole genes: {multi}/{total} = {100*multi/total:.1f}%"
              f"   ⟹ {multi}/{len(dropped)} = {100*multi/len(dropped):.1f}% of ALL evictions")
        print(f"  distinct culprit loci {len(culprit)}; top-10 share "
              f"{100*sum(c for _, c in culprit.most_common(10))/total:.1f}% (diffuse, not a few monsters)")


if __name__ == '__main__':
    main()
