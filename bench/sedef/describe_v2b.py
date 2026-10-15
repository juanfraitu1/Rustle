#!/usr/bin/env python3
"""Descriptive (no bar): which Soto sets make up V2b's expansion_not_family in each arm, and how many of them equal a
`<family>.unique` bucket. Reuses bench/o1_hier/hier_match.py unchanged (same inputs as the registered V2b command).

    v2b_detail.py <arm dir>"""
import collections
import importlib.util
import sys

spec = importlib.util.spec_from_file_location(
    "hier_match", "/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto/bench/o1_hier/hier_match.py")
hm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hm)

arm = sys.argv[1]
ext = "/mnt/linuxdisk/tmp/hier_flag/soto_external.tsv"
nodes, family_of = hm.load_nodes(f"{arm}/ours.hierarchy.members.tsv")
members = sorted(family_of)
sets, per_record, cov = hm.truth_sets(hm.load_truth(ext), members, hm.load_fold(f"{arm}/ours.loci.tsv"), None, "key")
ev = {k: v for k, v in sets.items() if len(v) >= 2}
print(arm, "mapped Soto sets with >= 2 members:", len(ev))

levels_of = collections.defaultdict(set)
names_of = collections.defaultdict(list)
for (level, nid), s in nodes.items():
    levels_of[s].add(level)
    names_of[s].append((level, nid))

hits = {}
for label, s in ev.items():
    ls = levels_of.get(s, set())
    if "expansion" in ls and "family" not in ls:
        hits[label] = sorted(nid for (lv, nid) in names_of[s] if lv == "expansion")

uniq = {k: v for k, v in hits.items() if all(n.endswith(".unique") for n in v)}
numbered = {k: v for k, v in hits.items() if not all(n.endswith(".unique") for n in v)}
print("expansion_not_family:", len(hits), " of which equal a .unique bucket:", len(uniq), " numbered expansions:", len(numbered))
for label in sorted(hits):
    print("  ", label, hits[label], "size", len(ev[label]))
