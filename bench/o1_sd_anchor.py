#!/usr/bin/env python3
"""Does SEGMENTAL-DUPLICATION membership separate O1's false merges from its true pairs?

WHY THIS IS THE RIGHT SHAPE. The 2026-08-20 synthesis: every statistic computed FROM THE TWO NODES
re-encodes node length and collapses when length-normalised (coverage on all four substrate x
denominator cells, block count, junction crossing, read tiling). Only a statistic referencing
something OUTSIDE the pair escapes. An SD catalog is exactly that — and unlike the gene annotation it
is not entailed by a gene-based truth predicate.

THE TEST. For each catalog pair (a, b), is there an SD pair whose two units land one on a's locus and
one on b's? A duplicated gene family should sit inside a segmental duplication. A repeat-bridged false
merge should not: the bridge is a mobile element, not a duplicated genomic segment.

⚠ NOT independent evidence. An SD is DEFINED by high identity over >=1 kb, so SD membership is itself
a homology statement. What it adds is GENOMIC CONTIGUITY and the duplication UNIT's boundaries —
things the exon-sum substrate cannot see. Read it as a better substrate, not an oracle.

Substrate: SEDEF final.bed on mGorGor1 (253,030 pairs) + the 494-family catalog the frozen arms use.
UNIT = PAIR. GGO. FP arm 14 (NOT held out), TP arm 150 (load-bearing). T8: offline.
"""
import csv, collections, sys

SD  = "/mnt/c/Users/jfris/Desktop/final.bed"
D   = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
A   = "/home/juanfra/winloci_scratch/o1_antifp"
MIN_ID  = float(sys.argv[1]) if len(sys.argv) > 1 else 0.90   # canonical SD identity floor
MIN_LEN = 1000                                                # canonical SD length floor


def load_sd():
    """chrom -> sorted list of (start, end, partner_chrom, p_start, p_end, identity), both directions"""
    by = collections.defaultdict(list)
    kept = tot = 0
    for line in open(SD):
        f = line.rstrip("\n").split("\t")
        if len(f) < 34:
            continue
        tot += 1
        try:
            ident = float(f[-1])
        except ValueError:
            continue
        a0, a1, b0, b1 = int(f[1]), int(f[2]), int(f[4]), int(f[5])
        if ident < MIN_ID or (a1 - a0) < MIN_LEN or (b1 - b0) < MIN_LEN:
            continue
        kept += 1
        by[f[0]].append((a0, a1, f[3], b0, b1, ident))
        by[f[3]].append((b0, b1, f[0], a0, a1, ident))
    for c in by:
        by[c].sort()
    return by, kept, tot


def main():
    sd, kept, tot = load_sd()
    print(f"SD catalog: {tot} rows -> {kept} pass identity >= {MIN_ID} and both units >= {MIN_LEN} bp", flush=True)

    cp = {f"{r['family_id']}~{r['copy_idx']}": (r["chrom"], int(r["start"]), int(r["end"]))
          for r in csv.DictReader(open(f"{D}/GGO_gwcat.copies.tsv"), delimiter="\t")}
    fp = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{A}/fp14_detail.tsv"), delimiter="\t")]
    tp = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{A}/tp150_detail.tsv"), delimiter="\t")]

    def linked(a, b):
        """is there an SD pair placing one unit on a and the other on b? returns best identity or None"""
        ca, sa, ea = cp[a]; cb, sb, eb = cp[b]
        best = None
        for u0, u1, pc, p0, p1, ident in sd.get(ca, []):
            if u0 >= ea:
                break
            if u1 <= sa:
                continue
            if pc == cb and p0 < eb and sb < p1:
                best = ident if best is None else max(best, ident)
        return best

    out = {}
    for lab, pairs in (("FP", fp), ("TP", tp)):
        hits = 0; idents = []
        for a, b in pairs:
            if a not in cp or b not in cp:
                continue
            v = linked(a, b)
            if v is not None:
                hits += 1; idents.append(v)
        out[lab] = (hits, len(pairs), idents)

    print("\n=== IS THE PAIR SPANNED BY A SEGMENTAL DUPLICATION? (unit = pair) ===")
    print(f"  {'arm':<6} {'SD-linked':>12} {'rate':>8}   median SD identity")
    for lab in ("FP", "TP"):
        h, n, idn = out[lab]
        med = sorted(idn)[len(idn)//2] if idn else float("nan")
        print(f"  {lab:<6} {h:>5}/{n:<6} {h/n:>8.4f}   {med:.4f}" if idn
              else f"  {lab:<6} {h:>5}/{n:<6} {h/n:>8.4f}   n/a")
    hf, nf, _ = out["FP"]; ht, nt, _ = out["TP"]
    print(f"\n  As an ADMISSION criterion (require SD support): keeps {ht}/{nt} = {ht/nt:.4f} of true pairs,")
    print(f"  admits {hf}/{nf} = {hf/nf:.4f} of false merges.")
    print(f"  As a VETO (reject when SD-supported): would reject {hf} FP and cost {ht} TP — nonsense unless FP >> TP.")


if __name__ == "__main__":
    main()
