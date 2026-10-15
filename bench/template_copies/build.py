#!/usr/bin/env python3
"""The copy's loci from its own reads (spec section 1, step 5).

The projection chooses the region and the reads; the reads give every coordinate. A read's exon blocks that overlap the
projected exons are kept; overlapping kept blocks of different reads merge into one observed exon; two consecutive blocks
of one read (no dropped block between them) link their observed exons. Each connected component is one locus: its
transcript is the path when the component is a simple path (links exactly between neighbouring observed exons), else
each distinct read chain is a transcript and chains contained in another are dropped. A component touched by fewer than
`min_reads` reads is not emitted (the read-support rule applied per locus; plan ruling R1). Nothing is filled from the
template and no boundary tolerance is used.
"""
import bisect
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "figures"))
import _liftoff as LF  # noqa: E402
from templates import merge  # noqa: E402


def _runs(blocks_, proj):
    """A read's kept blocks as runs of consecutive kept blocks: [[block, ...], ...]."""
    runs, prev = [], None
    for i, b in enumerate(blocks_):
        if LF.iv_inter([b], proj) > 0:
            if prev is not None and prev == i - 1:
                runs[-1].append(b)
            else:
                runs.append([b])
            prev = i
    return runs


def _contained(a, b):
    """True when tuple a is a contiguous sub-sequence of tuple b (a != b)."""
    return len(a) < len(b) and any(b[i:i + len(a)] == a for i in range(len(b) - len(a) + 1))


def build(reads, proj, min_reads=2):
    """reads = [(name, exon blocks)], proj = projected exons. -> loci, sorted by start:
    [{'transcripts': [[(s, e), ...], ...], 'reads': int, 'exons': [(s, e), ...]}]."""
    runs = [(name, run) for name, bl in reads for run in _runs(bl, proj)]
    obs = merge([b for _, run in runs for b in run])
    starts = [s for s, _ in obs]
    node = lambda b: bisect.bisect_right(starts, b[0]) - 1  # noqa: E731  (the observed exon containing block b)
    parent = list(range(len(obs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    edges, chains = set(), []
    for name, run in runs:
        chain = []
        for b in run:
            n = node(b)
            if chain and chain[-1] != n:
                edges.add((chain[-1], n))
                parent[find(n)] = find(chain[-1])
            if not chain or chain[-1] != n:
                chain.append(n)
        chains.append((name, tuple(chain)))
    comp_nodes, comp_reads, comp_chains = (collections.defaultdict(set) for _ in range(3))
    for n in range(len(obs)):
        comp_nodes[find(n)].add(n)
    for name, chain in chains:
        c = find(chain[0])
        comp_reads[c].add(name)
        comp_chains[c].add(chain)
    out = []
    for c, nodes in comp_nodes.items():
        if len(comp_reads[c]) < min_reads:
            continue
        ordered = sorted(nodes)
        cedges = {e for e in edges if find(e[0]) == c}
        if cedges == {(ordered[i], ordered[i + 1]) for i in range(len(ordered) - 1)}:
            txs = [ordered]
        else:
            cs = sorted(comp_chains[c])
            txs = [list(ch) for ch in cs if not any(_contained(ch, o) for o in cs)]
        out.append({"transcripts": [[obs[n] for n in t] for t in txs], "reads": len(comp_reads[c]),
                    "exons": [obs[n] for n in ordered]})
    return sorted(out, key=lambda lc: lc["exons"][0])
