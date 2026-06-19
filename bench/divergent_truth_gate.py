#!/usr/bin/env python3
"""Gate for the divergent-paralog track: how much INDEPENDENT paralog truth actually covers de-novo locus
pairs, broken out by core coverage (the divergent population is the LOW-core tail).

Truth tiers, all independent of our RNA k-mer feature:
  - gene-name ROOT (RefSeq curated names; same root = candidate paralog) -- LARGE but noisy for mega-families.
  - LOC genes carry NO name truth (unique numeric ids) -- the coverage wall.
Reads only existing artifacts (no network, no big compute). Output = the go/no-go numbers.
"""
import csv, re, bisect, collections, sys

ANNOT = "/home/juanfra/winloci_scratch/annot_intervals.tsv"
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
EDGES = "denovo_family_edges.tsv"

def gene_root(name):
    """Crude paralog root: strip a trailing copy-number (digits + optional single letter). APOBEC3D->APOBEC,
    RFPL2->RFPL, ZNF131->ZNF. LOC genes have no root (unique ids) -> return None (no name truth)."""
    if name.startswith("LOC"):
        return None
    m = re.match(r"^([A-Za-z][A-Za-z-]*?)\d+[A-Z]?$", name)
    return m.group(1) if m else name  # genes with no numeric suffix keep their full symbol

# --- genes per chrom, sorted by start ---
genes = collections.defaultdict(list)
for row in csv.DictReader(open(ANNOT), delimiter="\t"):
    genes[row["chrom"]].append((int(row["start"]), int(row["end"]), row["gene"], row["biotype"]))
for c in genes:
    genes[c].sort()
starts = {c: [g[0] for g in v] for c, v in genes.items()}

# root -> how many distinct genes carry it (to flag mega-family "blob" roots that are not clean paralogs)
root_size = collections.Counter()
for c in genes:
    for _, _, name, _ in genes[c]:
        r = gene_root(name)
        if r:
            root_size[r] += 1

def map_locus(chrom, s, e):
    """gene with the largest overlap of [s,e); None if no overlap."""
    if chrom not in genes:
        return None
    arr, st = genes[chrom], starts[chrom]
    i = bisect.bisect_right(st, e) - 1
    best, bov = None, 0
    while i >= 0 and arr[i][1] >= s - 200_000:  # walk back over genes that could overlap
        gs, ge, name, bt = arr[i]
        ov = min(e, ge) - max(s, gs)
        if ov > bov:
            bov, best = ov, (name, bt)
        if gs < s - 200_000:
            break
        i -= 1
    return best if bov > 0 else None

# --- locus coords from meta ---
coord = {}
for row in csv.DictReader(open(META), delimiter="\t"):
    coord[row["id"]] = (row["chrom"], int(row["start"]), int(row["end"]))

# --- classify every accepted edge pair, bucketed by core coverage ---
MEGA = 10  # a root carried by > MEGA genes is a blob (ZNF/OR...), not clean paralog truth
buckets = [(0.13, 0.15), (0.15, 0.20), (0.20, 0.30), (0.30, 1.01)]
def bkt(x):
    for lo, hi in buckets:
        if lo <= x < hi:
            return (lo, hi)
    return (0.13, 0.15)

stat = collections.defaultdict(lambda: collections.Counter())
n_loci_named = n_loci_loc = n_loci_unmapped = 0
seen = set()
for row in csv.DictReader(open(EDGES), delimiter="\t"):
    cr = float(row["core_recip"])
    b = bkt(cr)
    g = {}
    ok = True
    for k in ("a", "b"):
        lid = row[k]
        if lid not in coord:
            ok = False
            break
        g[k] = map_locus(*coord[lid])
        if lid not in seen:
            seen.add(lid)
            if g[k] is None:
                n_loci_unmapped += 1
            elif g[k][0].startswith("LOC"):
                n_loci_loc += 1
            else:
                n_loci_named += 1
    if not ok:
        continue
    ga, gb = g["a"], g["b"]
    stat[b]["pairs"] += 1
    if ga is None or gb is None:
        stat[b]["unmapped"] += 1
        continue
    ra, rb = gene_root(ga[0]), gene_root(gb[0])
    if ra is None or rb is None:
        stat[b]["has_LOC"] += 1
    elif ra == rb:
        stat[b]["same_root"] += 1
        if root_size[ra] <= MEGA:
            stat[b]["same_root_CLEAN"] += 1  # clean Tier-2 truth: small named family, same root
        else:
            stat[b]["same_root_blob"] += 1
    else:
        stat[b]["diff_named"] += 1  # both named, different root -> needs Compara to label

print(f"=== locus->gene coverage (distinct loci in edges) ===")
tot = n_loci_named + n_loci_loc + n_loci_unmapped
print(f"named={n_loci_named} ({100*n_loci_named/tot:.0f}%)  LOC={n_loci_loc} ({100*n_loci_loc/tot:.0f}%)  unmapped={n_loci_unmapped} ({100*n_loci_unmapped/tot:.0f}%)")
print(f"\n=== accepted edge pairs by core_recip bucket (the divergent population is the LOW buckets) ===")
print(f"{'bucket':<14}{'pairs':>7}{'same_root_CLEAN':>17}{'same_root_blob':>16}{'diff_named':>12}{'has_LOC':>9}{'unmapped':>10}")
for b in buckets:
    s = stat[b]
    print(f"{str(b):<14}{s['pairs']:>7}{s['same_root_CLEAN']:>17}{s['same_root_blob']:>16}{s['diff_named']:>12}{s['has_LOC']:>9}{s['unmapped']:>10}")
tot_clean = sum(stat[b]["same_root_CLEAN"] for b in buckets)
low_clean = stat[(0.13,0.15)]["same_root_CLEAN"] + stat[(0.15,0.20)]["same_root_CLEAN"]
print(f"\nCLEAN gene-name-truth pairs total: {tot_clean}  | in the LOW-core (<0.20, most divergent) bins: {low_clean}")
