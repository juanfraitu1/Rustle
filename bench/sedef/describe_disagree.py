#!/usr/bin/env python3
"""Descriptive breakdown (no bar) of where the sd-scan duplicon cover differs from the DupMasker cover, per half.

    v2_disagree.py <dupmasker dir> <sdscan dir>

Reuses bench/paf_duplicons/score_units.py's covers / family_score / half unchanged."""
import collections
import importlib.util
import sys

spec = importlib.util.spec_from_file_location(
    "score_units", "/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto/bench/paf_duplicons/score_units.py")
su = importlib.util.module_from_spec(spec)
spec.loader.exec_module(su)

truth, ours = sys.argv[1:3]
hd, rd = su._members(f"{truth}/ours.hierarchy.members.tsv")
hp, rp = su._members(f"{ours}/ours.hierarchy.members.tsv")
key = lambda h, r: (r[h["family"]], r[h["member"]], r[h["expansions"]])  # noqa: E731
dc, pc = su.covers([key(hd, r) for r in rd]), su.covers([key(hp, r) for r in rp])
assert sorted(dc) == sorted(pc)
print("families (>= 2 members):", len(dc), " NPIP family", su.NPIP_FAMILY, "is in half:", su.half(su.NPIP_FAMILY))

for which in ("heldout", "dev"):
    fams = sorted(f for f in dc if su.half(f) == which)
    cat = collections.Counter()
    members_by_cat = collections.Counter()
    sc_by_cat = collections.defaultdict(list)
    for f in fams:
        d_units = [k for k in dc[f] if k != "unique"]
        p_units = [k for k in pc[f] if k != "unique"]
        s = su.family_score(dc[f], pc[f])
        n = len(set().union(*dc[f].values()))
        if s == 1.0:
            c = "identical"
        elif d_units and not p_units:
            c = "DupMasker has expansions, sd-scan all unique"
        elif p_units and not d_units:
            c = "sd-scan has expansions, DupMasker all unique"
        elif d_units and p_units:
            c = "both have expansions, differ"
        else:
            c = "neither has expansions, differ"
        cat[c] += 1
        members_by_cat[c] += n
        sc_by_cat[c].append(s)
    print(f"\n== {which}: {len(fams)} families")
    for c, v in sorted(cat.items(), key=lambda x: -x[1]):
        print(f"  {v:4d}  members {members_by_cat[c]:5d}  mean score {sum(sc_by_cat[c]) / v:.4f}  {c}")
    nd = sum(1 for f in fams if len([k for k in dc[f] if k != 'unique']) > 0)
    npp = sum(1 for f in fams if len([k for k in pc[f] if k != 'unique']) > 0)
    print(f"  families with >= 1 non-unique block: DupMasker {nd}, sd-scan {npp}")
    nzero = sum(1 for f in fams if su.family_score(dc[f], pc[f]) == 0.0)
    print(f"  families scoring exactly 0.0: {nzero}")
