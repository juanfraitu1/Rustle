#!/usr/bin/env python3
"""Diff graph edges (adjacency) RU vs ST, keyed by bundle envelope.
Converts the edge-parity verdict from 'proved by reduction' to 'measured'.
Usage: python3 bench/edge_diff.py /tmp/ru_edge.jsonl /tmp/st_edge.jsonl"""
import json, sys
from collections import defaultdict

def load(path):
    edges = {}   # (start,end,strand) -> set of "FROM>TO"
    nodes = {}   # (start,end,strand) -> frozenset of "A-B" node tokens
    with open(path) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except ValueError:
                continue
            step = r.get("step")
            if step not in ("graph_edge", "graphnode_list"):
                continue
            key = (r["start"], r["end"], r["strand"])
            if step == "graph_edge":
                e = r.get("payload", {}).get("edges", "")
                edges[key] = set(x for x in e.split(",") if x)
            else:
                nstr = r.get("payload", {}).get("nodes", "")
                toks = set(m.split(":")[0] for m in nstr.split(",") if m.strip())
                nodes[key] = frozenset(toks)
    return edges, nodes

def edge_class(e):
    if ">" not in e:
        return "malformed"
    f, t = e.split(">", 1)
    if f == "SRC" or t == "SNK" or f == "SNK" or t == "SRC":
        return "source_sink"
    try:
        fa, fb = map(int, f.split("-"))
        ta, tb = map(int, t.split("-"))
    except ValueError:
        return "malformed"
    return "junction" if ta > fb + 1 else "collinear"

def main():
    ru_e, ru_n = load(sys.argv[1] if len(sys.argv) > 1 else "/tmp/ru_edge.jsonl")
    st_e, st_n = load(sys.argv[2] if len(sys.argv) > 2 else "/tmp/st_edge.jsonl")
    shared = set(ru_e) & set(st_e)
    identical_node = [k for k in shared if ru_n.get(k) == st_n.get(k)]
    print(f"edge events: RU {len(ru_e)} / ST {len(st_e)} / shared {len(shared)}")
    print(f"identical-node-set shared bundles: {len(identical_node)} (others have node diffs)")

    idn_edge_id = sum(1 for k in identical_node if ru_e[k] == st_e[k])
    print(f"  of identical-node bundles, IDENTICAL edge sets: {idn_edge_id}/{len(identical_node)}")
    all_id = sum(1 for k in shared if ru_e[k] == st_e[k])
    print(f"all shared bundles with IDENTICAL edge sets: {all_id}/{len(shared)}")

    st_only_class = defaultdict(int)   # ST has, RU lacks (RU-MISSING)
    ru_only_class = defaultdict(int)   # RU has, ST lacks (RU-EXTRA)
    ru_missing_junc_examples = []
    ru_missing_junc_idnode = 0
    for k in shared:
        for e in (st_e[k] - ru_e[k]):
            c = edge_class(e); st_only_class[c] += 1
            if c == "junction":
                if k in set(identical_node):
                    ru_missing_junc_idnode += 1
                if len(ru_missing_junc_examples) < 15:
                    ru_missing_junc_examples.append((k in set(identical_node), k, e))
        for e in (ru_e[k] - st_e[k]):
            ru_only_class[edge_class(e)] += 1

    print(f"RU-MISSING edges (ST has, RU lacks) by class: {dict(st_only_class)}")
    print(f"  >>> RU MISSING junction-spanning edges total: {st_only_class.get('junction',0)} "
          f"(of which on identical-NODE bundles: {ru_missing_junc_idnode}) — reduction predicts ~0")
    print(f"RU-EXTRA edges (RU has, ST lacks) by class: {dict(ru_only_class)}")
    print("examples of RU-MISSING junction edges (idnode?, bundle, edge):")
    for idn, k, e in ru_missing_junc_examples:
        if edge_class(e) == "junction":
            print(f"    idnode={idn} {k} {e}")

if __name__ == "__main__":
    main()
