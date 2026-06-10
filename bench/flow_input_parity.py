#!/usr/bin/env python3
"""Layered flow-input parity harness.

Partitions rustle's canonical-only intron chains (chains canonical extracts that StringTie does
not) into GRAPH-ROOTED (the local graph differs from ST → architectural ceiling, not fixable by
the flow algorithm) vs FLOW-ROOTED (rustle has all of ST's junctions in the region → the divergence
is purely flow path-selection → the flow-port's addressable target).

Method: for each canonical-only chain, take its bundle region (the rustle graph bundle containing
it). If StringTie has good-junctions in that region that rustle LACKS, the local graph differs
→ graph-rooted. If rustle has every ST junction in the region, the graphs match → flow-rooted.

Usage:
  python3 bench/flow_input_parity.py <ru_canon.gtf> <st.gtf> <ru_struct.jsonl> <st_struct.jsonl>
"""
import json, re, sys, collections

def safe(p):
    for l in open(p):
        try:
            yield json.loads(l)
        except ValueError:
            continue

def good_juncs(path):
    """Set of (start, end) good-junction coords (strand emit differs across tools → coords only)."""
    s = set()
    for o in safe(path):
        if o.get('step') == 'good_junction':
            s.add((o['start'], o['end']))
    return s

def bundles(path):
    """List of (start, end) bundle spans from graphnode_list events."""
    out = []
    for o in safe(path):
        if o.get('step') == 'graphnode_list':
            out.append((o['start'], o['end']))
    return out

def chains(path):
    """Return {(strand, intron_chain): sorted_exon_list}."""
    tx = collections.defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = {}
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2:
            continue
        ch = tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))
        out[(strand[t], ch)] = ex
    return out

def main():
    ru_gtf, st_gtf, ru_struct, st_struct = sys.argv[1:5]
    ru_c = chains(ru_gtf); st_c = chains(st_gtf)
    canon_only = {k: v for k, v in ru_c.items() if k not in st_c}
    ru_j = good_juncs(ru_struct); st_j = good_juncs(st_struct)

    print(f"flow-input parity: rustle good-juncs={len(ru_j)}  ST good-juncs={len(st_j)}  "
          f"(rustle subset of ST: {ru_j <= st_j})")
    print(f"canonical-only chains to partition: {len(canon_only)}\n")

    # Precise criterion: a chain is GRAPH-ROOTED if ST has a good-junction (that rustle LACKS)
    # whose start OR end lands strictly INSIDE one of the chain's exons — i.e. ST's graph would
    # split that exon, making this exact chain unbuildable in ST's graph. Such chains are
    # artifacts of rustle's coarser (missing-junction) graph, not flow path-selection.
    # FLOW-ROOTED: all the chain's exons survive intact in ST's graph (no internal ST junction
    # rustle lacks) → ST could build this chain too but its flow selected differently → the
    # flow-port's addressable target.
    st_only_j = st_j - ru_j  # ST junctions rustle lacks
    st_only_coords = sorted({c for jj in st_only_j for c in jj})
    import bisect

    def splits_an_exon(exons):
        """count ST-only junction endpoints landing strictly inside any exon."""
        n = 0
        for (es, ee) in exons:
            i = bisect.bisect_right(st_only_coords, es)
            while i < len(st_only_coords) and st_only_coords[i] < ee:
                n += 1
                i += 1
        return n

    graph_rooted = 0; flow_rooted = 0
    splits_total = 0
    details = []
    for (strand, chain), exons in canon_only.items():
        nsplit = splits_an_exon(exons)
        if nsplit > 0:
            graph_rooted += 1
            splits_total += nsplit
            details.append((strand, exons[0][0], exons[-1][1], len(exons), nsplit))
        else:
            flow_rooted += 1

    print(f"=== PARTITION of the {len(canon_only)} canonical-only chains "
          f"(precise: ST-only junction lands inside the chain's exon) ===")
    print(f"  GRAPH-ROOTED (an ST junction rustle lacks lands inside one of the chain's exons): "
          f"{graph_rooted}  ({100*graph_rooted/max(1,len(canon_only)):.0f}%)")
    print(f"     -> rustle's & ST's graphs differ in node structure here (different flow INPUTS);")
    print(f"        rooted in the junction-acceptance divergence, NOT the flow algorithm -> NOT flow-port-fixable")
    print(f"  FLOW-ROOTED  (chain's exons intact in ST's graph; same local node structure): "
          f"{flow_rooted}  ({100*flow_rooted/max(1,len(canon_only)):.0f}%)")
    print(f"     -> same graph inputs; ST's flow selected differently -> the flow-port's addressable target")
    if graph_rooted:
        details.sort(key=lambda d: -d[4])
        print(f"\n  graph-rooted chains have {splits_total} exon-splitting ST junctions total")
        print("  worst (strand span nexons exon_splits):")
        for s, lo, hi, nx, sp in details[:8]:
            print(f"    {s} {lo}-{hi}  nexons={nx} exon_splits={sp}")

if __name__ == '__main__':
    main()
