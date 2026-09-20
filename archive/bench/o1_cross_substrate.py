#!/usr/bin/env python3
"""Does O1's definition describe the GENOME, or the sample it was run on?

The shipped 627-family catalog was built from GGO_ds.bam = OR6737 TESTIS -- a different animal from the
assembly (KB3781). This compares it against a depth-matched catalog from KB3781 FIBROBLAST (1,627,629 vs
1,628,629 primaries), same reference, same binary, same parameters. Exactly two things differ: ANIMAL and
TISSUE.

⚠ COUNTS ARE NOT THE READOUT. Testis transcribes far more of the genome, so fewer fibroblast loci is
expected biology. Judging a definition on how many families each sample yields would repeat the error
this project has made three times (judging a change to what a NODE is on node-level counts). The
question is whether the two agree ON THE LOCI THEY SHARE -- the "same loci" discipline that showed DNA
and RNA give an identical partition 7/7.

⚠ AND AGREEMENT NEEDS ITS COMPARATOR. Two partitions of the same loci agree somewhat by chance, and
this project killed four headlines in one week to rates quoted without one. The null here permutes copy
labels WITHIN the shared locus set while preserving each catalog's family SIZE distribution (a
count-matched null proves nothing; the size distribution is what must match).
"""
import bisect, collections, csv, random, sys

T = "/mnt/linuxdisk/home/juanfraitu/o1_reps2"          # testis, determinism-gated == shipped
F = "/mnt/linuxdisk/home/juanfraitu/o1_replicate"      # fibroblast, depth-matched
TG = "/mnt/linuxdisk/home/juanfraitu/o1_gw/ggo_gw.copies.tsv"


def nodes(p):
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(p), delimiter="\t"):
        by[r["chrom"]].append((int(r["start"]), int(r["end"]), int(r["idx"]), int(r["degree"])))
    for c in by:
        by[c].sort()
    return by


def match(a, b, recip=0.50):
    """Reciprocal-overlap match between two rep sets. Reciprocal, not any-overlap: a 200 kb rep
    swallowing a 2 kb one is not the same locus."""
    out = {}
    for c, va in a.items():
        vb = b.get(c)
        if not vb:
            continue
        starts = [x[0] for x in vb]
        for s, e, i, d in va:
            j0 = bisect.bisect_left(starts, e)
            best = None
            for k in range(max(0, j0 - 60), j0):
                s2, e2, i2, d2 = vb[k]
                ov = min(e, e2) - max(s, s2)
                if ov <= 0:
                    continue
                if ov >= recip * (e - s) and ov >= recip * (e2 - s2):
                    if best is None or ov > best[0]:
                        best = (ov, i2)
            if best:
                out[i] = best[1]
    return out


def edges(p):
    E = set()
    for r in csv.DictReader(open(p), delimiter="\t"):
        i, j = int(r["rep_i"]), int(r["rep_j"])
        E.add((min(i, j), max(i, j)))
    return E


def ari(lab_a, lab_b):
    """Adjusted Rand Index over items labelled in both."""
    keys = [k for k in lab_a if k in lab_b]
    n = len(keys)
    if n < 2:
        return float("nan")
    tab = collections.Counter((lab_a[k], lab_b[k]) for k in keys)
    ra = collections.Counter(lab_a[k] for k in keys)
    rb = collections.Counter(lab_b[k] for k in keys)
    c2 = lambda x: x * (x - 1) / 2
    sij = sum(c2(v) for v in tab.values())
    sa = sum(c2(v) for v in ra.values()); sb = sum(c2(v) for v in rb.values())
    exp = sa * sb / c2(n); mx = (sa + sb) / 2
    return (sij - exp) / (mx - exp) if mx != exp else float("nan")


def main():
    ta, fa = nodes(f"{T}/dump/ggo.nodes.tsv"), nodes(f"{F}/dump/fibro.nodes.tsv")
    nt = sum(len(v) for v in ta.values()); nf = sum(len(v) for v in fa.values())
    m = match(ta, fa)
    print(f"reps: testis {nt}   fibroblast {nf}")
    print(f"  reciprocal-50% matched loci: {len(m)}  "
          f"= {len(m)/nt:.4f} of testis, {len(m)/nf:.4f} of fibroblast\n")

    Et, Ef = edges(f"{T}/dump/ggo.edges.tsv"), edges(f"{F}/dump/fibro.edges.tsv")
    inv = {v: k for k, v in m.items()}
    Et_s = {(a, b) for a, b in Et if a in m and b in m}
    Ef_s = {(min(inv[a], inv[b]), max(inv[a], inv[b])) for a, b in Ef if a in inv and b in inv}
    inter = Et_s & Ef_s; union = Et_s | Ef_s
    print(f"E_r edges among SHARED loci: testis {len(Et_s)}   fibroblast {len(Ef_s)}")
    print(f"  intersection {len(inter)}   union {len(union)}   "
          f"Jaccard {len(inter)/len(union) if union else float('nan'):.4f}")
    if Et_s:
        print(f"  of testis edges on shared loci, recovered by fibroblast: "
              f"{len(inter)}/{len(Et_s)} = {len(inter)/len(Et_s):.4f}")

    def labels(path, byc):
        lab = {}
        for r in csv.DictReader(open(path), delimiter="\t"):
            c, s, e = r["chrom"], int(r["start"]), int(r["end"])
            for (s2, e2, i, d) in byc.get(c, ()):
                if s2 == s and e2 == e:
                    lab[i] = r["family_id"]; break
        return lab
    lt = labels(TG, ta); lf = labels(f"{F}/fibro_gwcat.copies.tsv", fa)
    lt_s = {k: v for k, v in lt.items() if k in m}
    lf_s = {inv[k]: v for k, v in lf.items() if k in inv}
    both = set(lt_s) & set(lf_s)
    print(f"\ncatalog COPIES on shared loci: testis {len(lt_s)}  fibroblast {len(lf_s)}  both {len(both)}")
    if len(both) >= 2:
        a = {k: lt_s[k] for k in both}; b = {k: lf_s[k] for k in both}
        obs = ari(a, b)
        print(f"  ARI of the induced partitions: {obs:.4f}")
        rng = random.Random(101); null = []
        for _ in range(200):
            vals = list(b.values()); rng.shuffle(vals)
            null.append(ari(a, dict(zip(b.keys(), vals))))
        null.sort(); mu = sum(null) / len(null)
        ge = sum(1 for x in null if x >= obs)
        print(f"  size-preserving label-permutation null: mean {mu:.4f}, "
              f"95th pct {null[int(.95*len(null))]:.4f}, p = {(ge+1)/(len(null)+1):.4f}")


if __name__ == "__main__":
    main()
