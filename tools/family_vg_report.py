#!/usr/bin/env python3
"""
family_vg_report.py — build a "family variation graph" from annotated
transcripts of a gene family.

Usage:
    python3 family_vg_report.py <GFF> <pattern> [--tol=500] [--min-introns=5]

Takes gene/mRNA/exon records from a GFF and finds genes whose description
field matches <pattern> (e.g. 'golgin subfamily A member 6'). Clusters them
by Needleman-Wunsch chain alignment into sub-families. Emits the consensus
intron-length chain per sub-family and the per-member deviations.

This is analogous to what a tool like vg giraffe does conceptually — build
a graph representing the variation across a related set of paths — but
operates at the intron-length level instead of the nucleotide level.
"""
import re
import sys
import argparse
from collections import defaultdict


def nw_align(a, b, bp_penalty=0.1, gap_penalty=500, bp_cap=2000):
    """Needleman-Wunsch on intron-length sequences."""
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return -(max(n, m) * gap_penalty), []
    dp = [[0.0] * (m + 1) for _ in range(n + 1)]
    bt = [[None] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = -i * gap_penalty; bt[i][0] = 'U'
    for j in range(1, m + 1):
        dp[0][j] = -j * gap_penalty; bt[0][j] = 'L'
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pen = min(abs(a[i-1]-b[j-1]) * bp_penalty, bp_cap)
            d = dp[i-1][j-1] - pen
            u = dp[i-1][j] - gap_penalty
            l = dp[i][j-1] - gap_penalty
            best = max(d, u, l)
            dp[i][j] = best
            bt[i][j] = 'D' if best == d else ('U' if best == u else 'L')
    matches = []
    i, j = n, m
    while i > 0 or j > 0:
        if bt[i][j] == 'D':
            matches.append((i-1, j-1)); i -= 1; j -= 1
        elif bt[i][j] == 'U':
            i -= 1
        else:
            j -= 1
    matches.reverse()
    return dp[n][m], matches


def load_family(gff_path, pattern):
    genes = []
    tx_parent = {}
    exons_by_tx = defaultdict(list)
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fs = line.rstrip('\n').split('\t')
            if len(fs) < 9: continue
            if fs[2] == 'gene':
                m_name = re.search(r'Name=([^;]*)', fs[8])
                m_id = re.search(r'ID=([^;]*)', fs[8])
                m_desc = re.search(r'description=([^;]*)', fs[8])
                if not m_desc or pattern not in m_desc.group(1): continue
                genes.append((fs[0], int(fs[3]), int(fs[4]), fs[6],
                              m_id.group(1), m_name.group(1) if m_name else ''))
            elif fs[2] == 'mRNA':
                tx_parent[re.search(r'ID=([^;]*)', fs[8]).group(1)] = \
                    re.search(r'Parent=([^;]*)', fs[8]).group(1)
            elif fs[2] == 'exon':
                exons_by_tx[re.search(r'Parent=([^;]*)', fs[8]).group(1)].append(
                    (int(fs[3]), int(fs[4])))
    chains = {}
    for chrom, s, e, strand, gid, name in genes:
        best = None
        for tx, par in tx_parent.items():
            if par != gid: continue
            ex = sorted(exons_by_tx[tx])
            if len(ex) < 2: continue
            ch = [ex[i+1][0]-ex[i][1]-1 for i in range(len(ex)-1)]
            if best is None or len(ch) > len(best):
                best = ch
        if best:
            if strand == '-': best = list(reversed(best))
            chains[name] = (chrom, s, strand, best)
    return chains


def cluster(chains, score_threshold, min_introns):
    names = sorted(chains.keys())
    names = [n for n in names if len(chains[n][3]) >= min_introns]
    parent = {n: n for n in names}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, a in enumerate(names):
        for b in names[i+1:]:
            s, _ = nw_align(chains[a][3], chains[b][3])
            if s >= score_threshold:
                ra, rb = find(a), find(b)
                if ra != rb: parent[ra] = rb
    groups = defaultdict(list)
    for n in names: groups[find(n)].append(n)
    return list(groups.values())


def consensus(chains, members):
    """Align all members to a reference (the one with most introns),
    then emit the consensus chain with per-position variance."""
    ref = max(members, key=lambda m: len(chains[m][3]))
    ref_chain = chains[ref][3]
    # Map each aligned position of each member onto the reference
    per_pos = [[] for _ in range(len(ref_chain))]
    for m in members:
        if m == ref:
            for i, v in enumerate(ref_chain):
                per_pos[i].append(v)
        else:
            _, matches = nw_align(ref_chain, chains[m][3])
            for (ri, mi) in matches:
                per_pos[ri].append(chains[m][3][mi])
    return ref, ref_chain, per_pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gff')
    ap.add_argument('pattern', help="gene description substring (e.g. 'golgin subfamily A member 6')")
    ap.add_argument('--score', type=float, default=-500,
                    help='alignment score threshold for clustering (default: -500, more negative = looser)')
    ap.add_argument('--min-introns', type=int, default=5)
    args = ap.parse_args()

    chains = load_family(args.gff, args.pattern)
    if not chains:
        print(f"No genes matched pattern: {args.pattern}", file=sys.stderr)
        sys.exit(1)
    print(f"Loaded {len(chains)} family members from annotation")
    print()

    groups = cluster(chains, args.score, args.min_introns)
    multi = [g for g in groups if len(g) > 1]
    singletons = [g[0] for g in groups if len(g) == 1]
    print(f"=== Sub-families: {len(multi)} (with {sum(len(g) for g in multi)} members) + "
          f"{len(singletons)} singletons ===")
    print()

    for i, grp in enumerate(sorted(multi, key=lambda g: -len(g))):
        print(f"SUB-FAMILY {i+1}: {len(grp)} members")
        for m in sorted(grp, key=lambda n: (chains[n][0], chains[n][1])):
            c, s, st, ch = chains[m]
            print(f"    {m:<22s} {c:<14s} {s:>10d} {st}   {len(ch)} introns")
        ref, ref_chain, per_pos = consensus(chains, grp)
        conserved = sum(1 for p in per_pos if len(p) == len(grp) and max(p) - min(p) <= 20)
        print(f"    Consensus (relative to {ref}, {len(ref_chain)} introns, "
              f"{conserved} strictly conserved within 20bp):")
        for pos, vals in enumerate(per_pos):
            if not vals: continue
            mn, mx = min(vals), max(vals)
            mean = sum(vals) / len(vals)
            tag = "CONSERVED" if mx - mn <= 20 else "VARIABLE" if mx - mn <= 200 else "DIVERGED"
            spread = f"{mn:5d}-{mx:5d}" if mn != mx else f"{mn:5d}"
            print(f"      i{pos+1:02d}  n={len(vals)}/{len(grp)}  mean={int(mean):5d}  range={spread:>13s}  {tag}")
        print()

    if singletons:
        print(f"\n=== Singletons (no aligned partner at score threshold {args.score}) ===")
        for m in sorted(singletons, key=lambda n: (chains[n][0], chains[n][1])):
            c, s, st, ch = chains[m]
            print(f"    {m:<22s} {c:<14s} {s:>10d} {st}   {len(ch)} introns")


if __name__ == '__main__':
    main()
