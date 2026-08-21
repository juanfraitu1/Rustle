#!/usr/bin/env python3
"""Adjudicate the 114 questionable γ splits.

WHY THESE 114. Measured 2026-08-21: of 1,138 SD-supported node pairs that should be one family, γ
loses 246 = 21.6% — 3.0x what the coverage clause loses. But 132 of those lie inside the 530-copy
hairball and are FORCED cuts, where γ had to partition a blob no definition could emit whole. The
remaining 114 are outside it: γ chose to split them with room to do otherwise. They are the largest
unexamined object in O1 and they sit in the bigger of the two recoverable loss channels.

THE QUESTION. Is γ over-splitting real families, or correctly separating DIFFERENT genes that happen
to be co-duplicated inside one segmental duplication? An SD can carry two unrelated genes, so
SD support alone cannot tell them apart.

ADJUDICATION. Map each locus to its best-overlapping annotated gene and compare:

  same gene symbol          -> one gene's copies were split           OVER-SPLIT
  (a LOC id counts for EQUALITY but never for the stem test — LOC numbering is positional)
  same symbol STEM          -> paralogues (GSTM1/GSTM2)               OVER-SPLIT
  different unrelated genes -> co-duplicated distinct genes           γ CORRECT
  unnamed / LOC only        -> cannot adjudicate                      UNRESOLVED

⚠ The annotation is legitimate here: it produced neither the γ partition nor the SD calls, so this is
evaluation against evidence the thing being tested did not generate. ⚠ The stem heuristic (strip
trailing digits) is crude — GSTM1/GSTM2 share a stem, but so would unrelated genes numbered alike.
Every OVER-SPLIT call is printed so it can be eyeballed rather than trusted.
"""
import bisect, collections, csv, re, sys

G   = "/mnt/linuxdisk/home/juanfraitu/o1_gw"
D   = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
SD  = "/mnt/c/Users/jfris/Desktop/final.bed"
GFF = f"{D}/GGO_genomic.gff"


def main():
    fam, loc = {}, {}
    for r in csv.DictReader(open(f"{G}/ggo_gw.copies.tsv"), delimiter="\t"):
        k = f"{r['family_id']}~{r['copy_idx']}"
        fam[k] = r["family_id"]; loc[k] = (r["chrom"], int(r["start"]), int(r["end"]))

    E = set()
    for paf, mid in ((f"{G}/gw_ava_asm20.paf", 0.80), (f"{G}/gw_ava_sens.paf", 0.60)):
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[0] == f[5] or f[0] not in fam or f[5] not in fam:
                continue
            ql, qs, qe, tl = float(f[1]), float(f[2]), float(f[3]), float(f[6])
            de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
            i = (1.0 - de) if de is not None else int(f[9]) / max(int(f[10]), 1)
            if i >= mid and (qe - qs) / max(min(ql, tl), 1.0) >= 0.50:
                E.add((min(f[0], f[5]), max(f[0], f[5])))

    par = {}
    def find(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in E:
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb
    comp = collections.defaultdict(list)
    for k in fam:
        comp[find(k)].append(k)
    hairball = set(max(comp.values(), key=len))

    by = collections.defaultdict(list)
    for k, (c, s, e) in loc.items():
        by[c].append((s, e, k))
    for c in by:
        by[c].sort()
    st = {c: [x[0] for x in v] for c, v in by.items()}
    def hits(c, s, e):
        v = by.get(c)
        if not v:
            return []
        i = bisect.bisect_left(st[c], e); out = []
        for j in range(max(0, i - 500), i):
            s0, e0, k = v[j]
            if e0 > s:
                out.append(k)
        return out

    sdp = set()
    for line in open(SD):
        f = line.rstrip("\n").split("\t")
        if len(f) < 34:
            continue
        try:
            idn = float(f[-1])
        except ValueError:
            continue
        a0, a1, b0, b1 = int(f[1]), int(f[2]), int(f[4]), int(f[5])
        if idn < 0.90 or a1 - a0 < 1000 or b1 - b0 < 1000:
            continue
        for x in hits(f[0], a0, a1):
            for y in hits(f[3], b0, b1):
                if x != y and fam[x] != fam[y]:
                    sdp.add((min(x, y), max(x, y)))

    gsplit = [p for p in sdp if p in E]
    outside = [p for p in gsplit if not (p[0] in hairball and p[1] in hairball)]
    print(f"SD-supported cross-family pairs passing E_r: {len(gsplit)}")
    print(f"  inside the 530-copy hairball (forced):     {len(gsplit)-len(outside)}")
    print(f"  OUTSIDE it — the questionable set:         {len(outside)}\n", flush=True)

    genes = collections.defaultdict(list)
    for line in open(GFF):
        if line[0] == "#":
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.search(r"(?:^|;)gene=([^;]+)", f[8])
        if m:
            genes[f[0]].append((int(f[3]) - 1, int(f[4]), m.group(1)))
    for c in genes:
        genes[c].sort()
    def sym(k):
        c, s, e = loc[k]
        best, bov = None, 0
        for g0, g1, n in genes.get(c, []):
            if g0 >= e:
                break
            ov = min(e, g1) - max(s, g0)
            if ov > bov:
                best, bov = n, ov
        return best
    stem = lambda x: re.sub(r"\d+$", "", x) if x else None

    cls = collections.Counter(); over = []
    for a, b in outside:
        sa, sb = sym(a), sym(b)
        # ⚠ A LOC identifier is still an IDENTITY, so `sa == sb` adjudicates it. What a LOC id cannot
        # support is the STEM heuristic — LOC numbering is positional, not family-based, so LOC123456
        # and LOC123457 are neighbours, not paralogues. Excluding LOC entirely (the first version of
        # this script) threw away 97.4% of the adjudicable set.
        both_loc = (sa or "").startswith("LOC") and (sb or "").startswith("LOC")
        if sa is None or sb is None:
            cls["UNRESOLVED (no gene)"] += 1
        elif sa == sb:
            cls["OVER-SPLIT (same gene)"] += 1; over.append((a, b, sa, sb))
        elif both_loc:
            # ⚠⚠ NOT evidence of correctness. RefSeq assigns a distinct LOC id PER LOCUS, so two
            # copies of one UNNAMED family necessarily carry different LOC ids. This class cannot
            # distinguish "two copies of one unnamed family" from "two different unnamed genes",
            # so it is UNINFORMATIVE. Counting it as γ-correct inflates the exoneration ~7x.
            cls["UNINFORMATIVE (different LOC ids)"] += 1
        elif stem(sa) == stem(sb) and len(stem(sa)) >= 3:
            cls["OVER-SPLIT (same stem)"] += 1; over.append((a, b, sa, sb))
        else:
            cls["γ CORRECT (different genes)"] += 1

    print("=== adjudication of the questionable γ splits ===")
    n = len(outside)
    for k, v in cls.most_common():
        print(f"  {k:<32} {v:>4}/{n}  = {v/n:.4f}")
    adj = sum(v for k, v in cls.items() if k.startswith("OVER-SPLIT") or k.startswith("γ CORRECT"))
    if adj:
        o = cls["OVER-SPLIT (same gene)"] + cls["OVER-SPLIT (same stem)"]
        print(f"\n  among ADJUDICABLE pairs ({adj}): over-split {o}/{adj} = {o/adj:.4f}")
    print("\n  every OVER-SPLIT call, for eyeballing:")
    for a, b, sa, sb in over[:25]:
        print(f"    {fam[a]:<10} {sa:<14} | {fam[b]:<10} {sb}")
    if len(over) > 25:
        print(f"    ... {len(over)-25} more")


if __name__ == "__main__":
    main()
