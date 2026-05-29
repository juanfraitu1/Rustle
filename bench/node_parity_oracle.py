#!/usr/bin/env python3
"""Layer-2.5 ORACLE: does graph-node parity predict transfrag-chain parity?

Partitions bundles into (node-set matches ST) vs (doesn't), then measures the
transfrag_pre_depl Rustle-only rate in each partition. If node-matching bundles have
~0 Rustle-only transfrags, node-parity is the lever -> build Layer 2.5 first.
Inputs: /tmp/ru_both.jsonl /tmp/st_both.jsonl (graphnode_list + transfrag_pre_depl).
"""
import json, re, collections

def load(path):
    gn = {}                              # bundle key -> set of node intervals
    tf = collections.defaultdict(set)    # bundle key -> set of (strand,intron-chain)
    rows = []
    for line in open(path):
        try: rows.append(json.loads(line))
        except Exception: continue
    for e in rows:
        if e.get("step") != "graphnode_list": continue
        p = e.get("payload", e)
        iv = set()
        for tok in p.get("nodes", "").split(","):
            m = re.match(r"(\d+)-(\d+)", tok.strip())
            if m: iv.add((int(m.group(1)), int(m.group(2))))
        gn[(e.get("start"), e.get("end"), e.get("strand"))] = iv
    for e in rows:
        if e.get("step") != "transfrag_pre_depl": continue
        p = e.get("payload", e)
        ich = p.get("introns", "").rstrip(",")
        if not ich: continue
        ts, te, tstr = e.get("start"), e.get("end"), e.get("strand")
        placed = False
        for bkey in gn:
            bs, be, bstr = bkey
            if bs is not None and ts is not None and bs <= ts and te <= be and tstr == bstr:
                tf[bkey].add((tstr, ich)); placed = True; break
        if not placed:
            tf[("orphan",)].add((tstr, ich))
    return gn, tf

ru_gn, ru_tf = load("/tmp/ru_both.jsonl")
st_gn, st_tf = load("/tmp/st_both.jsonl")
shared = set(ru_gn) & set(st_gn)
match_b    = [k for k in shared if ru_gn[k] == st_gn[k]]
mismatch_b = [k for k in shared if ru_gn[k] != st_gn[k]]

def ru_only_chains(bkeys):
    return sum(len(ru_tf.get(k, set()) - st_tf.get(k, set())) for k in bkeys)

total_ru_only = sum(len(ru_tf.get(k, set()) - st_tf.get(k, set())) for k in ru_tf)
print(f"shared bundles: {len(shared)}  node-MATCH: {len(match_b)}  node-MISMATCH: {len(mismatch_b)}")
print(f"Rustle-only transfrag chains, node-MATCH bundles:    {ru_only_chains(match_b)}")
print(f"Rustle-only transfrag chains, node-MISMATCH bundles: {ru_only_chains(mismatch_b)}")
print(f"Rustle-only transfrag chains, orphan/unplaced:       {len(ru_tf.get(('orphan',), set()))}")
print(f"total Rustle-only (all partitions):                  {total_ru_only}")
