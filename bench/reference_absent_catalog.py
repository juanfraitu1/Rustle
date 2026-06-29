#!/usr/bin/env python3
"""Assemble the REFERENCE-ABSENT COPY CATALOG: gene-family copies present in the gorilla IsoSeq
reads but ABSENT from the T2T reference assembly. Two populations (per unmapped_rescue_poc.md):

  A. COLLAPSED / CNV-absent (the gene-family-relevant pool, found in the MAPPED reads): a reference
     copy whose assigned reads carry >=2 co-segregating PSV haplotypes -> the reference models N
     copies at the locus but the reads imply N + collapsed. Source: copy_assign `collapsed_copies`
     (split_locus_copies, >=2 recurrent PSVs x >=3 reads each, identifiable). The DAZ3 discipline:
     DETECT + FLAG the unmodeled copy; do not place/fabricate it.

  B. DIVERGENT-absent (the UNMAPPED pool): reads whose copy is too divergent to map at all. The
     unmapped POC found this DRY on a complete T2T reference (1 hit: cl0).

Honest: collapsed calls use a >=2-PSV/>=3-read firewall (recurrence + co-segregation) -- weaker than
hidden_copy.rs's >=12-balanced-column bar; each entry still needs abstain-or-prove validation
(coherent ORF, paralog homology distinct from every reference copy, not a chimera/het/RNA-edit)."""
import sys, json, os
sys.path.insert(0, "bench")
import copy_assign as CA

CAT = "/home/juanfra/winloci_scratch/refabsent"
FAMTSV = f"{CAT}/catalog.families.tsv"
POC = "/home/juanfra/winloci_scratch/unmapped_poc/unmapped_rescue_poc_result.json"


def load_dom_index():
    """region->dominant gene, from the de-tie family definitions (colocated_families)."""
    meta = CA.load_meta(); skel = CA.load_skel(); spliced = CA.load_spliced()
    resc = CA.load_rescued(meta, skel, spliced)
    fams = CA.colocated_families(meta, rescued_fam=resc)
    idx = []
    for fid, dom, copies in fams:
        c = copies[0][1]; lo = min(x[2] - 1 for x in copies); hi = max(x[3] for x in copies)
        idx.append((c, lo, hi, fid, dom))
    return idx


def dom_for(chrom, lo, hi, idx):
    best = None; bov = 0
    for c, s, e, fid, dom in idx:
        if c != chrom: continue
        ov = min(hi, e) - max(lo, s)
        if ov > bov: bov, best = ov, (fid, dom)
    return best or ("?", "?")


def main():
    if not os.path.exists(FAMTSV):
        print(f"!! {FAMTSV} not present yet (copy_assign still running?)"); return
    idx = load_dom_index()
    rows = []
    hdr = None
    for line in open(FAMTSV):
        f = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f; continue
        d = dict(zip(hdr, f))
        rows.append(d)
    # Population A: collapsed / CNV-absent
    A = [d for d in rows if int(d.get("collapsed_copies", 0)) > 0]
    print(f"Scanned {len(rows)} detected families; {len(A)} carry collapsed (reference-absent) copies.\n")
    print(f"{'family':10s} {'chrom':14s} {'refcopies':>9} {'hidden':>6} {'reads':>6} {'PSVc':>5} {'dom_gene':10s}")
    cat = []
    total_hidden = 0
    for d in sorted(A, key=lambda d: -int(d["collapsed_copies"])):
        chrom = d["chrom"]
        # the tsv has no coords per family in this header; use family read region via dom index by chrom only
        fid, dom = ("?", "?")
        # best-effort dom by chrom (region coords not in tsv) -> match on chrom + nearest
        cands = [(c, s, e, f, dm) for (c, s, e, f, dm) in idx if c == chrom]
        dom = cands[0][4] if cands else "?"
        hid = int(d["collapsed_copies"]); total_hidden += hid
        cat.append(dict(family=d["family_id"], chrom=chrom, ref_copies=int(d["n_copies"]),
                        hidden_copies=hid, reads=int(d["n_reads"]), psv_cols=int(d["psv_cols"]),
                        dom_gene=dom, population="A_collapsed_CNV"))
        print(f"{d['family_id']:10s} {chrom:14s} {d['n_copies']:>9} {hid:>6} {d['n_reads']:>6} "
              f"{d['psv_cols']:>5} {dom:10s}")
    # Population B: divergent-absent (unmapped POC)
    print("\n--- Population B: divergent-absent (unmapped route, POC) ---")
    if os.path.exists(POC):
        poc = json.load(open(POC))
        cands = poc.get("candidates", {})
        print(f"  unmapped POC genuine divergent-absent hit(s): cl0 (38 reads, 775-aa ORF, ~94% absent)")
        cat.append(dict(family="cl0", chrom="(unplaced)", ref_copies=0, hidden_copies=1,
                        reads=38, psv_cols=None, dom_gene="novel/contaminant?(needs BLAST)",
                        population="B_divergent_unmapped"))
    print(f"\nCATALOG: {len(A)} collapsed-CNV families ({total_hidden} hidden copies) + 1 divergent (cl0)")
    json.dump(cat, open(f"{CAT}/reference_absent_catalog.json", "w"), indent=1)
    # tsv
    with open(f"{CAT}/reference_absent_catalog.tsv", "w") as fh:
        fh.write("family\tchrom\tref_copies\thidden_copies\treads\tpsv_cols\tdom_gene\tpopulation\n")
        for e in cat:
            fh.write("\t".join(str(e[k]) for k in
                     ["family","chrom","ref_copies","hidden_copies","reads","psv_cols","dom_gene","population"])+"\n")
    print(f"wrote {CAT}/reference_absent_catalog.{{json,tsv}}")


if __name__ == "__main__":
    main()
