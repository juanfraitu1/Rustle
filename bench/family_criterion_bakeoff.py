#!/usr/bin/env python3
"""Conflict-edge criterion bake-off: AS-tie vs NM-tie vs de-tie on a labelled family panel.

The advisor's concern: minimap2's AS is a composite of length + the aligner's gap/clip scoring model — it
should not be taken at face value. This measures whether a tie criterion built on raw-er, length-aware
signals (NM = edit distance, de = gap-compressed per-base divergence) reproduces or improves on AS-tie for
the conflict-graph family definition, on Compara-labelled ground truth.

For every read with alignments in two member spans we attribute each ALIGNMENT RECORD to the single
best-overlapping member (the faithful analog of the Rust `build_read_placements` best-overlap fix that
killed the nested-gene false positive), then test the read's two placements for a TIE under each criterion:

  AS-tie : min(AS)        >= AS_TIE * max(AS),  AS > 0            (the current Rust default, as_tie=0.9)
  NM-tie : |r_a - r_b| <= D_NM  AND  max(r_a, r_b) <= CEIL_NM     where r = NM / query_aligned_len
  de-tie : |de_a-de_b| <= D_DE  AND  max(de_a,de_b) <= CEIL_DE    de = minimap2 de:f tag

A pair is a conflict EDGE iff conflict-read count >= MIN_READS (matches the Rust min_reads=2 default).

Panel: 5 Compara paralog pairs + 7 domain-sharers (the labelled ground truth) + cross-chromosome paralogs
(separate-chrom; firing is a legitimate TP iff reads truly cross-map, silent is correct otherwise) +
co-located arrays (MAGEA/PRAMEF/XAGE; should form one connected conflict component).

Output: bench/family_criterion_bakeoff.tsv (per-pair, all criteria) + a per-criterion confusion summary.
"""
import collections
import json
import os
import sys

import pysam

BAM = "/home/juanfra/winloci_scratch/GGO.bam"
OUT = os.environ.get("BAKEOFF_OUT", "/mnt/c/Users/jfris/Desktop/Rustle/bench/family_criterion_bakeoff.tsv")

# tie parameters (env-overridable for the sensitivity sweep)
AS_TIE = float(os.environ.get("AS_TIE", 0.9))      # Rust default
D_NM = float(os.environ.get("D_NM", 0.005))        # NM-rate difference tolerance (0.5%)
CEIL_NM = float(os.environ.get("CEIL_NM", 0.05))   # a read "fits" a copy iff its NM-rate <= 5%
D_DE = float(os.environ.get("D_DE", 0.005))        # de difference tolerance
CEIL_DE = float(os.environ.get("CEIL_DE", 0.05))   # a read "fits" iff de <= 5%
MIN_READS = int(os.environ.get("MIN_READS", 2))    # Rust default min_reads

# gene -> (chrom, start, end), from GGO_genomic.gff
GENES = {
    "RFPL1": ("NC_086018.1", 27445867, 27472024), "RFPL2": ("NC_086018.1", 30208643, 30221673),
    "RFPL3": ("NC_086018.1", 30368560, 30376055), "APOBEC3D": ("NC_086018.1", 37023493, 37035944),
    "APOBEC3F": ("NC_086018.1", 37043278, 37058621), "RABL2B": ("NC_086018.1", 48818440, 48832011),
    "RABL2A": ("NC_073235.2", 15131653, 15147533), "PPARA": ("NC_086018.1", 44159462, 44252981),
    "CDPF1": ("NC_086018.1", 44250569, 44262082), "CASP8": ("NC_073235.2", 103850223, 103901051),
    "FLACC1": ("NC_073235.2", 103898904, 103977183), "CREB1": ("NC_073235.2", 110133100, 110220981),
    "METTL21A": ("NC_073235.2", 110198842, 110241136), "GPR39": ("NC_073235.2", 34964263, 35195860),
    "LYPD1": ("NC_073235.2", 35194087, 35224540), "GCA": ("NC_073235.2", 64785386, 64835227),
    "KCNH7": ("NC_073235.2", 64832333, 65309180), "ASDURF": ("NC_073235.2", 92106380, 92111318),
    "ASNSD1": ("NC_073235.2", 92111059, 92115776), "ZDHHC8": ("NC_086018.1", 17564727, 17581526),
    "CCDC188": ("NC_086018.1", 17579900, 17589084),
    # cross-chromosome paralogs (gene + LOC copy on a different chrom)
    "ADH5": ("NC_073227.2", 119368977, 119386893), "LOC101125574": ("NC_073229.2", 89958139, 89959453),
    "AK6": ("NC_073243.2", 40875017, 40892984), "LOC115934278": ("NC_073227.2", 47415194, 47415981),
    "ATP5MF": ("NC_073230.2", 89357879, 89371496), "LOC101136485": ("NC_073228.2", 27516532, 27516798),
    "CCDC196": ("NC_073239.2", 83376335, 83388628), "LOC129526440": ("NC_073238.2", 33024134, 33036450),
    "CNN2": ("NC_073244.2", 6807132, 6819584), "LOC129524764": ("NC_073237.2", 31580552, 31581609),
    "EEF1A1": ("NC_073229.2", 97608632, 97613071), "LOC109023808": ("NC_073243.2", 97380144, 97381766),
    "ELF2": ("NC_073227.2", 159850315, 159969202), "LOC115933780": ("NC_073235.2", 87996030, 87997236),
    # co-located array members
    "MAGEA9": ("NC_073247.2", 161458542, 161464324), "MAGEA4": ("NC_073247.2", 163590519, 163600784),
    "MAGEA10": ("NC_073247.2", 163809309, 163814522), "MAGEA12": ("NC_073247.2", 164405219, 164411360),
    "MAGEA1": ("NC_073247.2", 164861291, 164865952),
    "PRAMEF20": ("NC_073224.2", 226778521, 226789458), "PRAMEF17": ("NC_073224.2", 226805764, 226808744),
    "PRAMEF19": ("NC_073224.2", 226826959, 226830413), "PRAMEF12": ("NC_073224.2", 227814760, 227818916),
    "XAGE2": ("NC_073247.2", 62535784, 62542447), "XAGE5": ("NC_073247.2", 62936754, 62944018),
    "XAGE3": ("NC_073247.2", 62986053, 62991620),
}

# label, gene_a, gene_b — the labelled ground truth
PARALOG = [("APOBEC3D", "APOBEC3F"), ("RFPL1", "RFPL2"), ("RFPL1", "RFPL3"),
           ("RFPL2", "RFPL3"), ("RABL2A", "RABL2B")]
DOMAIN = [("ASDURF", "ASNSD1"), ("CASP8", "FLACC1"), ("CCDC188", "ZDHHC8"), ("CDPF1", "PPARA"),
          ("CREB1", "METTL21A"), ("GCA", "KCNH7"), ("GPR39", "LYPD1")]
CROSSCHROM = [("ADH5", "LOC101125574"), ("AK6", "LOC115934278"), ("ATP5MF", "LOC101136485"),
              ("CCDC196", "LOC129526440"), ("CNN2", "LOC129524764"), ("EEF1A1", "LOC109023808"),
              ("ELF2", "LOC115933780")]
ARRAYS = {
    "MAGEA": ["MAGEA9", "MAGEA4", "MAGEA10", "MAGEA12", "MAGEA1"],
    "PRAMEF": ["PRAMEF20", "PRAMEF17", "PRAMEF19", "PRAMEF12"],
    "XAGE": ["XAGE2", "XAGE5", "XAGE3"],
}


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def records(bam, gene):
    """qname -> list of dicts (one per non-supplementary alignment record overlapping the gene span)."""
    c, s, e = GENES[gene]
    out = collections.defaultdict(list)
    for r in bam.fetch(c, s, e):
        if r.is_supplementary or r.is_unmapped:
            continue
        qlen = r.query_alignment_length or 0
        if qlen <= 0:
            continue
        try:
            as_ = r.get_tag("AS")
        except KeyError:
            as_ = 0
        try:
            nm = r.get_tag("NM")
        except KeyError:
            nm = 0
        try:
            de = float(r.get_tag("de"))
        except KeyError:
            de = nm / qlen  # fall back to raw NM-rate if de absent
        out[r.query_name].append(dict(
            chrom=c, rs=r.reference_start, re=r.reference_end, ov=overlap(r.reference_start, r.reference_end, s, e),
            as_=as_, nm_rate=nm / qlen, de=de, mapq=r.mapping_quality, sec=r.is_secondary))
    return out


def best_for_member(recs, gene):
    """The single record best-overlapping `gene`'s span (the best-overlap attribution)."""
    c, s, e = GENES[gene]
    cand = [x for x in recs if x["chrom"] == c and overlap(x["rs"], x["re"], s, e) > 0]
    if not cand:
        return None
    return max(cand, key=lambda x: x["ov"])


def tie_as(pa, pb):
    hi, lo = max(pa["as_"], pb["as_"]), min(pa["as_"], pb["as_"])
    return hi > 0 and lo >= AS_TIE * hi


def tie_nm(pa, pb, d=D_NM, ceil=CEIL_NM):
    return abs(pa["nm_rate"] - pb["nm_rate"]) <= d and max(pa["nm_rate"], pb["nm_rate"]) <= ceil


def tie_de(pa, pb, d=D_DE, ceil=CEIL_DE):
    return abs(pa["de"] - pb["de"]) <= d and max(pa["de"], pb["de"]) <= ceil


def conflict_counts(bam, ga, gb):
    """For each read placed (best-overlap) on BOTH ga and gb via DISTINCT records, count ties per criterion."""
    ra, rb = records(bam, ga), records(bam, gb)
    shared = set(ra) & set(rb)
    n_as = n_nm = n_de = 0
    disagree = []   # reads where criteria disagree (for the workflow's investigation)
    for q in shared:
        pa = best_for_member(ra[q], ga)
        pb = best_for_member(rb[q], gb)
        if pa is None or pb is None:
            continue
        # the two placements must be DISTINCT alignment records (different locus), not one record
        # straddling two nested annotations: require different ref coordinates or different chrom.
        if pa["chrom"] == pb["chrom"] and abs(pa["rs"] - pb["rs"]) < 200:
            continue
        t_as, t_nm, t_de = tie_as(pa, pb), tie_nm(pa, pb), tie_de(pa, pb)
        n_as += t_as
        n_nm += t_nm
        n_de += t_de
        if len({t_as, t_nm, t_de}) > 1:
            disagree.append(dict(q=q, as_tie=t_as, nm_tie=t_nm, de_tie=t_de,
                                 as_a=pa["as_"], as_b=pb["as_"], nmr_a=round(pa["nm_rate"], 4),
                                 nmr_b=round(pb["nm_rate"], 4), de_a=pa["de"], de_b=pb["de"]))
    return dict(n_shared=len(shared), n_as=n_as, n_nm=n_nm, n_de=n_de, disagree=disagree)


def main():
    bam = pysam.AlignmentFile(BAM, "rb")
    rows = []
    all_disagree = {}
    for label, pairs in [("paralog", PARALOG), ("domain-sharer", DOMAIN), ("crosschrom", CROSSCHROM)]:
        for ga, gb in pairs:
            cc = conflict_counts(bam, ga, gb)
            same_chrom = GENES[ga][0] == GENES[gb][0]
            edge = {k: (cc[f"n_{k}"] >= MIN_READS) for k in ("as", "nm", "de")}
            rows.append(dict(label=label, a=ga, b=gb, same_chrom=same_chrom, **cc, edge=edge))
            if cc["disagree"]:
                all_disagree[f"{ga}~{gb}"] = cc["disagree"]
            print(f"{label:13} {ga:13}~{gb:13} chrom_same={int(same_chrom)} shared={cc['n_shared']:4} "
                  f"| AS={cc['n_as']:3}{'E' if edge['as'] else ' '} "
                  f"NM={cc['n_nm']:3}{'E' if edge['nm'] else ' '} "
                  f"de={cc['n_de']:3}{'E' if edge['de'] else ' '}")

    # arrays: pairwise among members -> union-find component per criterion
    print("\n=== co-located arrays (member-pair conflict graph) ===")
    array_rows = []
    for fam, members in ARRAYS.items():
        for crit in ("as", "nm", "de"):
            parent = {m: m for m in members}

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x
            nedge = 0
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    cc = conflict_counts(bam, members[i], members[j])
                    if cc[f"n_{crit}"] >= MIN_READS:
                        nedge += 1
                        ri, rj = find(members[i]), find(members[j])
                        if ri != rj:
                            parent[ri] = rj
            comps = len({find(m) for m in members})
            array_rows.append(dict(family=fam, criterion=crit, n_members=len(members),
                                   n_edges=nedge, n_components=comps, connected=(comps == 1)))
            print(f"  {fam:8} [{crit}] members={len(members)} edges={nedge} components={comps} "
                  f"{'CONNECTED' if comps == 1 else 'SPLIT'}")

    with open(OUT, "w") as fh:
        fh.write("label\tgene_a\tgene_b\tsame_chrom\tn_shared\tn_as\tn_nm\tn_de\tedge_as\tedge_nm\tedge_de\n")
        for r in rows:
            fh.write(f"{r['label']}\t{r['a']}\t{r['b']}\t{int(r['same_chrom'])}\t{r['n_shared']}\t"
                     f"{r['n_as']}\t{r['n_nm']}\t{r['n_de']}\t{int(r['edge']['as'])}\t{int(r['edge']['nm'])}"
                     f"\t{int(r['edge']['de'])}\n")
        fh.write("\n# arrays\nfamily\tcriterion\tn_members\tn_edges\tn_components\tconnected\n")
        for r in array_rows:
            fh.write(f"{r['family']}\t{r['criterion']}\t{r['n_members']}\t{r['n_edges']}\t"
                     f"{r['n_components']}\t{int(r['connected'])}\n")
    with open(OUT.replace(".tsv", "_disagreements.json"), "w") as fh:
        json.dump(all_disagree, fh, indent=1)

    # confusion summary on the labelled pairs (paralog-with-conflict = should fire; domain-sharer = should not)
    print("\n=== confusion on labelled pairs (RFPL excluded from 'should-fire': resolvable by design) ===")
    should_fire = {("APOBEC3D", "APOBEC3F"), ("RABL2A", "RABL2B")}
    should_not = set(DOMAIN)
    for crit in ("as", "nm", "de"):
        tp = sum(1 for r in rows if (r["a"], r["b"]) in should_fire and r["edge"][crit])
        fn = len(should_fire) - tp
        fp = sum(1 for r in rows if (r["a"], r["b"]) in should_not and r["edge"][crit])
        tn = len(should_not) - fp
        print(f"  {crit.upper():3}: TP={tp}/{len(should_fire)} FN={fn} | TN={tn}/{len(should_not)} FP={fp}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    sys.exit(main())
