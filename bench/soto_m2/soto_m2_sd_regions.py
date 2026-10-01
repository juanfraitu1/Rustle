#!/usr/bin/env python3
"""SD98 regions and their duplicons, for the meeting page's "SD strips" (bench/soto_m2/page).

For every SD98 region (Soto's ≥ 98%-identity segmental duplications, merged, CHM13 v1.0) that overlaps at least one gene
of the page's gene list, records the genes inside it (indices into families.json's `genes`) and the DupMasker duplicon
segments inside it, clipped to the region. Duplicon colours are DupMasker's own (the itemRgb column), so the strips match
the colouring used in the SD literature.

    python3 bench/soto_m2/soto_m2_sd_regions.py --data families.json \
        --sd98 sd98_v1.bed --duplicons chm13.draft_v1.0_plus38Y_dupmasker_colors.bed --out sd_regions.json

Output: {"regions": [[chr, start, end, [gene idx...], [[start, end, duplicon idx, colour idx]...]], ...],
         "dups": [duplicon IDs], "cols": ["#rrggbb", ...]}; coordinates 0-based half-open, as in the BEDs."""
import argparse
import bisect
import collections
import json


def read_bed3(path):
    by = collections.defaultdict(list)
    for line in open(path):
        if line.strip() and not line.startswith(("#", "track", "browser")):
            f = line.split("\t")
            by[f[0]].append((int(f[1]), int(f[2])))
    for v in by.values():
        v.sort()
    return by


def overlapping(ivs, starts, a, b):
    """Indices of the sorted, non-overlapping intervals ivs that overlap [a, b)."""
    k = max(bisect.bisect_right(starts, a) - 1, 0)
    out = []
    while k < len(ivs) and ivs[k][0] < b:
        if ivs[k][1] > a:
            out.append(k)
        k += 1
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="families.json from soto_m2_families.py (its `genes` list)")
    ap.add_argument("--sd98", required=True, help="merged SD98 regions, BED, CHM13 v1.0")
    ap.add_argument("--duplicons", required=True, help="DupMasker duplicons with itemRgb, BED9, CHM13 v1.0")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    genes = json.load(open(a.data))["genes"]
    sd = read_bed3(a.sd98)
    for c, v in sd.items():
        assert all(v[k][1] <= v[k + 1][0] for k in range(len(v) - 1)), f"SD98 regions overlap on {c}"
    starts = {c: [x for x, _ in v] for c, v in sd.items()}
    held = collections.defaultdict(list)
    for i, g in enumerate(genes):
        for k in overlapping(sd.get(g["c"], []), starts.get(g["c"], []), g["s"], g["e"]):
            held[(g["c"], k)].append(i)
    unplaced = sum(1 for i, g in enumerate(genes)
                   if not overlapping(sd.get(g["c"], []), starts.get(g["c"], []), g["s"], g["e"]))

    dup_id, col_id, dups, cols = {}, {}, [], []
    segs = collections.defaultdict(list)
    for line in open(a.duplicons):
        f = line.rstrip("\n").split("\t")
        c, s, e = f[0], int(f[1]), int(f[2])
        for k in overlapping(sd.get(c, []), starts.get(c, []), s, e):
            if (c, k) not in held:
                continue
            r0, r1 = sd[c][k]
            if f[3] not in dup_id:
                dup_id[f[3]] = len(dups)
                dups.append(f[3])
            rgb = "#" + "".join(f"{int(x):02x}" for x in f[8].split(","))
            if rgb not in col_id:
                col_id[rgb] = len(cols)
                cols.append(rgb)
            segs[(c, k)].append([max(s, r0), min(e, r1), dup_id[f[3]], col_id[rgb]])

    regions = []
    for (c, k) in sorted(held, key=lambda ck: (ck[0], ck[1])):
        r0, r1 = sd[c][k]
        regions.append([c, r0, r1, sorted(held[(c, k)], key=lambda i: genes[i]["s"]), sorted(segs[(c, k)])])
    json.dump({"regions": regions, "dups": dups, "cols": cols}, open(a.out, "w"), separators=(",", ":"))
    print(f"{len(regions)} SD98 regions hold {sum(len(r[3]) for r in regions)} gene placements "
          f"({len(genes)} genes, {unplaced} outside every region); {sum(len(r[4]) for r in regions)} duplicon segments, "
          f"{len(dups)} duplicons, {len(cols)} colours")


if __name__ == "__main__":
    main()
