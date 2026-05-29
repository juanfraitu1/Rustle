#!/usr/bin/env python3
"""Layer-2.5 gate: compare graph node sets (graphnode_list) Rustle vs StringTie.

Joins by bundle key (chrom,start,end,strand). For each bundle present in BOTH tools,
compares the set of node (start,end) intervals (coverage ignored). Reports how many
shared bundles have IDENTICAL node sets and the aggregate node-interval divergence.
Layer 2.5 is done when shared-bundle node sets match (mismatch -> 0 / accounted).
Inputs: /tmp/ru_gn.jsonl /tmp/st_gn.jsonl (capture with PARITY_FILTER_STEPS=graphnode_list).
"""
import json, sys, re

def load(path):
    d = {}  # (start,end,strand) -> set of (node_start,node_end)
    for line in open(path):
        try: e = json.loads(line)
        except Exception: continue
        if e.get("step") != "graphnode_list": continue
        p = e.get("payload", e)
        nodes_s = p.get("nodes", "")
        key = (e.get("start"), e.get("end"), e.get("strand"))
        intervals = set()
        for tok in nodes_s.split(","):
            tok = tok.strip()
            if not tok: continue
            m = re.match(r"(\d+)-(\d+)", tok)
            if m: intervals.add((int(m.group(1)), int(m.group(2))))
        d[key] = intervals
    return d

ru = load(sys.argv[1] if len(sys.argv) > 1 else "/tmp/ru_gn.jsonl")
stt = load(sys.argv[2] if len(sys.argv) > 2 else "/tmp/st_gn.jsonl")
shared = set(ru) & set(stt)
print(f"Rustle bundles: {len(ru)}  ST bundles: {len(stt)}  shared: {len(shared)}")
identical = sum(1 for k in shared if ru[k] == stt[k])
ru_extra = sum(len(ru[k] - stt[k]) for k in shared)
st_extra = sum(len(stt[k] - ru[k]) for k in shared)
print(f"shared bundles with IDENTICAL node sets: {identical}/{len(shared)}")
print(f"node intervals Rustle-extra: {ru_extra}  ST-extra: {st_extra}")
print(f"\nLAYER-2.5 GATE (shared bundles with mismatched node sets): {len(shared) - identical} (target 0)")
