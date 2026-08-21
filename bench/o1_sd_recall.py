#!/usr/bin/env python3
"""What does SD-anchoring ADD to E_r?

SD membership is an ADMISSION certificate with FP 0/14 at every identity floor 0.99-0.80
(`o1_sd_anchor.py`). Because the false merges are already excluded by it, an SD-unit DENOMINATOR
cannot add discrimination on the FP arm — the FPs never get scored. So the value has to be RECALL:
does SD support identify pairs of catalog copies that belong together and that `E_r` did NOT join?

METHOD. Drive from the SD catalog, not from the pairs: for each SD pair, find catalog copies
overlapping unit A and unit B. Every (copyA, copyB) so linked is a duplication-supported pair. Split:

  SAME family      -> E_r already found it (concordance)
  DIFFERENT family -> a CANDIDATE missed edge

⚠ A cross-family SD link is a CANDIDATE, not a confirmed miss: one segmental duplication can carry
two genuinely different genes, so co-membership of an SD does not by itself imply one gene family.
Reported as candidates, and stratified by SD identity so the recent, high-confidence end is separable.

Substrate: SEDEF final.bed on mGorGor1 + the CURRENT 627-family catalog. UNIT = PAIR. T8: offline.
"""
import bisect, collections, csv, sys

SD  = "/mnt/c/Users/jfris/Desktop/final.bed"
CAT = "/mnt/linuxdisk/home/juanfraitu/o1_gw/ggo_gw.copies.tsv"
MIN_ID  = float(sys.argv[1]) if len(sys.argv) > 1 else 0.90
MIN_LEN = 1000


def main():
    # ---- catalog copies, indexed per contig ----
    by = collections.defaultdict(list)
    fam_of = {}
    for r in csv.DictReader(open(CAT), delimiter="\t"):
        k = f"{r['family_id']}~{r['copy_idx']}"
        by[r["chrom"]].append((int(r["start"]), int(r["end"]), k))
        fam_of[k] = r["family_id"]
    for c in by:
        by[c].sort()
    starts = {c: [x[0] for x in v] for c, v in by.items()}
    print(f"catalog: {len(fam_of)} copies in {len({v for v in fam_of.values()})} families", flush=True)

    def hits(chrom, s, e):
        v = by.get(chrom)
        if not v:
            return []
        i = bisect.bisect_left(starts[chrom], e)
        out = []
        for j in range(max(0, i - 400), i):
            s0, e0, k = v[j]
            if e0 > s:
                out.append(k)
        return out

    linked_same = set(); linked_diff = set(); ident_of = {}
    rows = 0
    for line in open(SD):
        f = line.rstrip("\n").split("\t")
        if len(f) < 34:
            continue
        try:
            ident = float(f[-1])
        except ValueError:
            continue
        a0, a1, b0, b1 = int(f[1]), int(f[2]), int(f[4]), int(f[5])
        if ident < MIN_ID or (a1 - a0) < MIN_LEN or (b1 - b0) < MIN_LEN:
            continue
        rows += 1
        ha, hb = hits(f[0], a0, a1), hits(f[3], b0, b1)
        if not ha or not hb:
            continue
        for x in ha:
            for y in hb:
                if x == y:
                    continue
                p = (min(x, y), max(x, y))
                (linked_same if fam_of[x] == fam_of[y] else linked_diff).add(p)
                if ident > ident_of.get(p, 0):
                    ident_of[p] = ident
    linked_diff -= linked_same
    print(f"SD pairs used (id >= {MIN_ID}, both units >= {MIN_LEN} bp): {rows}\n")

    print("=== duplication-supported catalog copy pairs ===")
    print(f"  SAME family      (E_r already joined them): {len(linked_same)}")
    print(f"  DIFFERENT family (CANDIDATE missed edges):  {len(linked_diff)}")
    tot = len(linked_same) + len(linked_diff)
    if tot:
        print(f"  ⟹ {len(linked_diff)/tot:.4f} of duplication-supported pairs are NOT in one family")

    print("\n=== candidates stratified by SD identity (the recent end is the confident end) ===")
    print(f"  {'SD identity':<14} {'same-family':>12} {'DIFFERENT':>11}")
    bands = [(0.99, 1.01), (0.98, 0.99), (0.95, 0.98), (0.90, 0.95)]
    for lo, hi in bands:
        s = sum(1 for p in linked_same if lo <= ident_of.get(p, 0) < hi)
        d = sum(1 for p in linked_diff if lo <= ident_of.get(p, 0) < hi)
        print(f"  {lo:.2f}-{hi:.2f}      {s:>12} {d:>11}")

    fams = collections.Counter()
    for x, y in linked_diff:
        fams[tuple(sorted((fam_of[x], fam_of[y])))] += 1
    print(f"\n  distinct FAMILY PAIRS implicated by candidates: {len(fams)}")
    print(f"  top: {fams.most_common(5)}")


if __name__ == "__main__":
    main()
