#!/usr/bin/env python
"""O3 genetic-core funnel — the REPRODUCIBLE replacement for the old '~77' headline (whose only backing,
bench/P1_P4_RESULTS.md, does not exist on disk; adversarial review #4). Builds the defensible genetic
allele-specific-junction core by stacking three filters on the committed call tables, with a per-call audit
TSV so the number is reproducible from a clean checkout.

Filter stack (each removes a distinct confound):
  1. TRANSVERSION anchor (high_confidence==1 in asj_calls_verified.tsv) — an A<->C/A<->T/C<->G/G<->T anchor base
     cannot be an A->I RNA-edit site, so allele-specificity is GENETIC not editing.   475 -> 120
  2. NON-LOC* (copy-confound) — LOC* loci are unannotated paralog families where a 'SNP' may distinguish two
     COPIES masquerading as two alleles (allele-vs-copy needs DNA); exclude them.       120 -> 76
  3. STRAND-BIAS clean (sor_pass==1 in asj_calls_strandbias.tsv; longcallR-style SOR) — removes mapping/strand
     artifacts. PSMD2 & DAXX (the old 'mechanistically airtight' flagships) FAIL here.   76 -> 54

The genetic core is the 54 strand-bias-clean transversion non-LOC* calls. (The old '~77' = 76 + 1 'clean LOC*'
added from a paralog-masquerade scan whose output is not committed; that +1 is retracted here, and the SOR
filter — the project's own — further tightens 76 -> 54.)

Run: /home/juanfra/miniforge3/bin/python bench/asj_genetic_core.py
"""
import csv
import os

HERE = os.path.dirname(__file__)


def load(name):
    return list(csv.DictReader(open(os.path.join(HERE, name)), delimiter="\t"))


def main():
    ver = load("asj_calls_verified.tsv")      # has high_confidence (transversion)
    sb = load("asj_calls_strandbias.tsv")      # has sor_pass
    assert len(ver) == len(sb), "call tables differ in length"
    # the two tables are the SAME 475 calls in the SAME order (verified by full-key match)
    key = lambda r: (r["gene"], r["chrom"], r["anchor"], r["alleleX"], r["alleleY"], r["donor"], r["acceptor"])
    assert all(key(a) == key(b) for a, b in zip(ver, sb)), "call tables are not row-aligned"

    rows = []
    for v, s in zip(ver, sb):
        transversion = v["high_confidence"] == "1"
        non_loc = not v["gene"].startswith("LOC")
        sor_pass = s["sor_pass"] == "1"
        in_core = transversion and non_loc and sor_pass
        rows.append({**v, "non_loc_copy_confound_clean": int(non_loc),
                     "sor_pass": int(sor_pass), "sor": s["sor"], "in_genetic_core": int(in_core)})

    n_all = len(rows)
    n_trans = sum(r["high_confidence"] == "1" for r in rows)
    n_nonloc = sum(r["high_confidence"] == "1" and r["non_loc_copy_confound_clean"] for r in rows)
    core = [r for r in rows if r["in_genetic_core"]]
    genes = sorted({r["gene"] for r in core})

    out = os.path.join(HERE, "asj_genetic_core.tsv")
    cols = ["gene", "chrom", "anchor", "alleleX", "alleleY", "donor", "acceptor", "dPSI", "q",
            "high_confidence", "non_loc_copy_confound_clean", "sor", "sor_pass", "in_genetic_core"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("== O3 genetic-core funnel (reproducible) ==")
    print(f"  {n_all} candidate ASJ calls")
    print(f"  -> {n_trans} TRANSVERSION (anchor cannot be A->I editing)")
    print(f"  -> {n_nonloc} NON-LOC* (copy/paralog-masquerade excluded)")
    print(f"  -> {len(core)} STRAND-BIAS clean (sor_pass) = the GENETIC CORE, in {len(genes)} genes")
    print(f"  core genes: {', '.join(genes)}")
    # show the retired flagships
    for g in ("PSMD2", "DAXX"):
        hit = next((r for r in rows if r["gene"] == g and r["high_confidence"] == "1"), None)
        if hit:
            print(f"  [retired flagship] {g}: transversion=1 sor={hit['sor']} sor_pass={hit['sor_pass']} "
                  f"-> {'IN' if hit['in_genetic_core'] else 'EXCLUDED from'} core")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
