#!/usr/bin/env python3
"""Sweep the family-grouping identity floor (--min-identity) and learn what predicts recoverability of the
74 coverage-independent floor members. The advisor's premise: this Soto set is the hardest-but-not-impossible
regime, so the failures should stratify into learnable patterns rather than one uniform wall.

Two axes per floor member:
  sibling_id  = best genomic identity to a DIFFERENT member (the grouping-relevant divergence; recomputed
                excluding self, since the earlier best_sibling_id was self-contaminated)
  unique_own  = uniquely-placed (MAPQ>=60) reads on the member (distinguishability, from the decomposition)

Recovery is checked against catalogs re-run at min-identity 0.98 (baseline), 0.95, 0.90 on the ideal-coverage
(top-up) BAM, so coverage and threshold are separated. Precision cost = copies that hit no annotated member.

Run: /home/juanfra/miniforge3/bin/python bench/soto/soto_threshold_sweep_eval.py
"""
import csv
from collections import defaultdict

D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
CATS = [("0.98", f"{D}/soto_topup.copies.tsv"),
        ("0.95", f"{D}/soto_topup_mid095.copies.tsv"),
        ("0.90", f"{D}/soto_topup_mid090.copies.tsv")]
DECOMP = "bench/soto/soto_floor_decomposition.tsv"
PAF = "bench/soto/a119b_member_pairs.paf"
BED = "bench/soto/80_fams.chr.bed"


def load_copies(path):
    idx = defaultdict(list)
    try:
        r = csv.reader(open(path), delimiter="\t"); next(r, None)
        for row in r:
            if len(row) >= 6:
                try:
                    idx[row[3]].append((int(row[4]), int(row[5])))
                except ValueError:
                    pass
    except FileNotFoundError:
        return None
    return idx


def hit(idx, chrom, s, e):
    return idx is not None and any(not (a > e or b < s) for (a, b) in idx.get(chrom, ()))


def sibling_ids():
    best = {}
    for line in open(PAF):
        f = line.split("\t")
        if len(f) < 11 or f[0] == f[5]:      # exclude self-alignment
            continue
        al = int(f[10])
        idp = int(f[9]) / al if al else 0
        if idp > best.get(f[0], 0):
            best[f[0]] = idp
    return best


def load_members_bed():
    idx = defaultdict(list)
    for line in open(BED):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[0].startswith("chr"):
            idx[p[0]].append((int(p[1]), int(p[2])))
    return idx


def band(x):
    if x is None:
        return "no-sibling"
    if x >= 0.999:
        return ">=99.9% (K=0)"
    if x >= 0.98:
        return "98-99.9%"
    if x >= 0.95:
        return "95-98%"
    if x >= 0.90:
        return "90-95%"
    return "<90%"


def main():
    cats = [(t, load_copies(p)) for t, p in CATS]
    if any(idx is None for _, idx in cats):
        print("waiting: not all catalogs present yet:", [t for (t, _), (_, i) in zip(CATS, cats) if i is None])
        return
    sib = sibling_ids()
    floor = list(csv.DictReader(open(DECOMP), delimiter="\t"))

    rows = []
    for r in floor:
        chrom, s, e = r["chrom"], int(r["start"]), int(r["end"])
        key = f"{chrom}:{s+1}-{e}"
        sid = sib.get(key)
        uniq = int(r["unique_own_reads"])
        rec = {t: ("Y" if hit(idx, chrom, s, e) else "N") for t, idx in cats}
        rows.append({"fam": r["family_id"], "gene": r["gene"], "chrom": chrom, "start": s, "end": e,
                     "sibling_id": sid, "unique_own": uniq, "band": band(sid),
                     "cause": r["cause"].split(" (")[0], **{f"rec_{t}": rec[t] for t, _ in cats}})

    # recovery by identity band
    print("\n=== recovery of the 74 floor members as the grouping floor is lowered ===")
    tot = {t: sum(1 for r in rows if r[f"rec_{t}"] == "Y") for t, _ in CATS}
    print(f"  detected@0.98={tot['0.98']}   @0.95={tot['0.95']}   @0.90={tot['0.90']}   (of {len(rows)})")
    print("\n  by sibling-identity band  (n : recovered@0.98 / @0.95 / @0.90 ; median unique reads):")
    bands = ["no-sibling", ">=99.9% (K=0)", "98-99.9%", "95-98%", "90-95%", "<90%"]
    byb = defaultdict(list)
    for r in rows:
        byb[r["band"]].append(r)
    for b in bands:
        g = byb.get(b, [])
        if not g:
            continue
        r98 = sum(1 for r in g if r["rec_0.98"] == "Y")
        r95 = sum(1 for r in g if r["rec_0.95"] == "Y")
        r90 = sum(1 for r in g if r["rec_0.90"] == "Y")
        us = sorted(r["unique_own"] for r in g)
        print(f"    {b:16s} n={len(g):2d} : {r98} / {r95} / {r90}   ; med_uniq={us[len(us)//2]}")

    # distinguishability of the newly-recovered vs still-missing (at 0.90)
    newly = [r for r in rows if r["rec_0.90"] == "Y" and r["rec_0.98"] == "N"]
    still = [r for r in rows if r["rec_0.90"] == "N"]
    def med_u(g):
        u = sorted(x["unique_own"] for x in g); return u[len(u)//2] if u else 0
    print(f"\n  newly recovered @0.90 (n={len(newly)}): median unique reads={med_u(newly)}, "
          f"median sibling_id={sorted(r['sibling_id'] for r in newly if r['sibling_id'] is not None)[len(newly)//2] if newly else 'NA':.4}" if newly else f"  newly recovered @0.90: 0")
    print(f"  still missing @0.90 (n={len(still)}): median unique reads={med_u(still)}; "
          f"causes: {dict((c, sum(1 for r in still if r['cause']==c)) for c in set(r['cause'] for r in still))}")

    # precision proxy: copies hitting no annotated member (potential over-merge / discovery)
    memidx = load_members_bed()
    print("\n=== precision cost: total copies + off-annotation copies per threshold ===")
    for t, idx in cats:
        allc = [(c, a, b) for c in idx for (a, b) in idx[c]]
        off = sum(1 for (c, a, b) in allc if not hit(memidx, c, a, b))
        print(f"    {t}: total_copies={len(allc):4d}  off-annotation={off:4d} ({100*off/max(len(allc),1):.0f}%)")

    with open("bench/soto/soto_threshold_sweep.tsv", "w", newline="") as f:
        cols = ["fam", "gene", "chrom", "start", "end", "sibling_id", "unique_own", "band", "cause",
                "rec_0.98", "rec_0.95", "rec_0.90"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print("\n  wrote bench/soto/soto_threshold_sweep.tsv")


if __name__ == "__main__":
    main()
