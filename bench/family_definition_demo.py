#!/usr/bin/env python3
"""FORMAL DEFINITION DEMONSTRATION — find multi-copy gene families from IsoSeq reads, with an honest FP/FN ledger.

THE DEFINITION (read-conflict graph, de-tie criterion):
  Let V = candidate loci. Build an undirected graph G=(V,E). For an ordered locus pair (i,j), a read CONFLICTS
  on (i,j) iff it has a best-overlap alignment on i AND a distinct best-overlap alignment on j whose
  gap-compressed divergences are TIED:  |de_i - de_j| <= DELTA (=0.005)  AND  max(de_i,de_j) <= DE_MAX (=0.05)
  (the read fits both loci comparably, at the HiFi error floor). Edge (i,j) in E iff >= MIN_READS (=3) reads
  conflict on it. A MULTI-COPY GENE FAMILY is a connected component of G with >= 2 loci.

  This defines the IDENTIFIABILITY-RELEVANT family: the set of loci among which reads are genuinely confused —
  exactly the unit the copy-assignment problem operates on. It is NOT a claim about evolutionary paralogy
  (RNA cannot recover that). A real paralog whose reads place uniquely (RFPL) is correctly NOT a family here.

PROPERTIES (what makes it a defensible formal object, not a tuned heuristic):
  (P1) Domain-sharer exclusion BY CONSTRUCTION: a read over a shared exon has ONE best-overlap placement, so it
       contributes no conflict pair -> nested/domain-sharing genes never form an edge. (No threshold involved.)
  (P2) Exact decomposition: reads never conflict across components, so the assignment problem separates exactly
       across families with no information lost (Canzar's conflict graph).
  (P3) One data-derived parameter: DELTA is the sequencing-divergence floor (HiFi de ~0.005-0.02); it is not a
       similarity threshold tuned per dataset. de is minimap2's raw divergence, not the aligner's composite AS.

This script instantiates the definition on a labelled panel of named gene families + decoys, recomputes the
conflict graph at the EXACT operating point from GGO.bam, reports the families found as components, and
tabulates every TP/TN/FP/FN. For contrast it also reports what the SIMILARITY definition (minimizer-Jaccard)
and the AS-tie criterion would do — to show our definition avoids their false positives.
"""
import collections
import json

import pysam

BAM = "/home/juanfra/winloci_scratch/GGO.bam"
OUT_TSV = "/mnt/c/Users/jfris/Desktop/Rustle/bench/family_definition_demo.tsv"
OUT_JSON = "/mnt/c/Users/jfris/Desktop/Rustle/bench/family_definition_demo.json"

# operating point (the shipped Rust de-tie defaults) — hardcoded so this demo cannot be clobbered by a sweep.
DELTA, DE_MAX, MIN_READS = 0.005, 0.05, 3
AS_TIE = 0.9  # legacy criterion, for the contrast column

GENES = {  # locus -> (chrom, start, end), from GGO_genomic.gff
    "RABL2A": ("NC_073235.2", 15131653, 15147533), "RABL2B": ("NC_086018.1", 48818440, 48832011),
    "APOBEC3D": ("NC_086018.1", 37023493, 37035944), "APOBEC3F": ("NC_086018.1", 37043278, 37058621),
    "RFPL1": ("NC_086018.1", 27445867, 27472024), "RFPL2": ("NC_086018.1", 30208643, 30221673),
    "RFPL3": ("NC_086018.1", 30368560, 30376055),
    "AK6": ("NC_073243.2", 40875017, 40892984), "LOC115934278": ("NC_073227.2", 47415194, 47415981),
    "CCDC196": ("NC_073239.2", 83376335, 83388628), "LOC129526440": ("NC_073238.2", 33024134, 33036450),
    "EEF1A1": ("NC_073229.2", 97608632, 97613071), "LOC109023808": ("NC_073243.2", 97380144, 97381766),
    "CNN2": ("NC_073244.2", 6807132, 6819584), "LOC129524764": ("NC_073237.2", 31580552, 31581609),
    "ASDURF": ("NC_073235.2", 92106380, 92111318), "ASNSD1": ("NC_073235.2", 92111059, 92115776),
    "CASP8": ("NC_073235.2", 103850223, 103901051), "FLACC1": ("NC_073235.2", 103898904, 103977183),
    "CREB1": ("NC_073235.2", 110133100, 110220981), "METTL21A": ("NC_073235.2", 110198842, 110241136),
    "GCA": ("NC_073235.2", 64785386, 64835227), "KCNH7": ("NC_073235.2", 64832333, 65309180),
    "GPR39": ("NC_073235.2", 34964263, 35195860), "LYPD1": ("NC_073235.2", 35194087, 35224540),
    "MAGEA9": ("NC_073247.2", 161458542, 161464324), "MAGEA4": ("NC_073247.2", 163590519, 163600784),
    "MAGEA10": ("NC_073247.2", 163809309, 163814522), "MAGEA1": ("NC_073247.2", 164861291, 164865952),
    # MAGEA region DE-NOVO expressed sub-loci (the vertices the shipped pipeline actually builds the conflict
    # graph over — coordinates from the detect_and_assign smoke run; the annotated MAGEA genes above are the
    # WRONG vertices for a tandem array because reads cross-map between these sub-loci, not the annotations).
    "MAGd0a": ("NC_073247.2", 161251228, 161257000), "MAGd0b": ("NC_073247.2", 161458538, 161464324),
    "MAGd1a": ("NC_073247.2", 161266095, 161270014), "MAGd1b": ("NC_073247.2", 161413981, 161453484),
    "MAGd2a": ("NC_073247.2", 164381222, 164384848), "MAGd2b": ("NC_073247.2", 164442447, 164446101),
    "MAGd3a": ("NC_073247.2", 164397061, 164401095), "MAGd3b": ("NC_073247.2", 164426194, 164430228),
}

# Ground-truth candidate families/decoys: (name, regime, [member loci], is_identifiability_family).
# is_identifiability_family = SHOULD our definition link them (reads genuinely cross-map)? — the truth label.
PANEL = [
    ("RABL2", "recent paralog, separate chrom", ["RABL2A", "RABL2B"], True),
    ("AK6", "high-identity cross-chrom copy", ["AK6", "LOC115934278"], True),
    ("CCDC196", "high-identity cross-chrom copy", ["CCDC196", "LOC129526440"], True),
    ("APOBEC3", "diverged paralog (resolvable)", ["APOBEC3D", "APOBEC3F"], False),
    ("RFPL", "paralog, reads place uniquely", ["RFPL1", "RFPL2", "RFPL3"], False),
    ("EEF1A1_retro", "processed-pseudogene retrocopy", ["EEF1A1", "LOC109023808"], False),
    ("CNN2_retro", "processed-pseudogene retrocopy", ["CNN2", "LOC129524764"], False),
    ("ASDURF_ds", "domain-sharer (nested genes)", ["ASDURF", "ASNSD1"], False),
    ("CASP8_ds", "domain-sharer", ["CASP8", "FLACC1"], False),
    ("CREB1_ds", "domain-sharer", ["CREB1", "METTL21A"], False),
    ("GCA_ds", "domain-sharer", ["GCA", "KCNH7"], False),
    ("GPR39_ds", "domain-sharer", ["GPR39", "LYPD1"], False),
    # co-located tandem array, DE-NOVO sub-loci vertices (the shipped pipeline's actual conflict families):
    ("MAGEA_dn0", "co-located array (de-novo loci)", ["MAGd0a", "MAGd0b"], True),
    ("MAGEA_dn1", "co-located array (de-novo loci)", ["MAGd1a", "MAGd1b"], True),
    ("MAGEA_dn2", "co-located array (de-novo loci)", ["MAGd2a", "MAGd2b"], True),
    ("MAGEA_dn3", "co-located array (de-novo loci)", ["MAGd3a", "MAGd3b"], True),
    # the ANNOTATED MAGEA genes are individually resolvable (reads do NOT cross-map between them) -> correctly
    # NOT one identifiability-family at the annotated-gene vertex set. (Contrast with the de-novo loci above.)
    ("MAGEA_annot", "annotated array genes (resolvable)", ["MAGEA9", "MAGEA4", "MAGEA10", "MAGEA1"], False),
]


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def records(bam, gene):
    c, s, e = GENES[gene]
    out = collections.defaultdict(list)
    for r in bam.fetch(c, s, e):
        if r.is_supplementary or r.is_unmapped:
            continue
        qlen = r.query_alignment_length or 0
        if qlen <= 0:
            continue
        out[r.query_name].append(dict(
            chrom=c, rs=r.reference_start, re=r.reference_end,
            ov=overlap(r.reference_start, r.reference_end, s, e),
            de=float(r.get_tag("de")) if r.has_tag("de") else r.get_tag("NM") / qlen,
            as_=r.get_tag("AS") if r.has_tag("AS") else 0, mapq=r.mapping_quality))
    return out


def best(recs, gene):
    c, s, e = GENES[gene]
    cand = [x for x in recs if x["chrom"] == c and overlap(x["rs"], x["re"], s, e) > 0]
    return max(cand, key=lambda x: x["ov"]) if cand else None


def pair_counts(bam, ga, gb):
    """conflict-read counts under de-tie and AS-tie + shared-read total, with best-overlap + distinct-record."""
    ra, rb = records(bam, ga), records(bam, gb)
    shared = set(ra) & set(rb)
    n_de = n_as = 0
    both_mapq0 = 0
    for q in shared:
        pa, pb = best(ra[q], ga), best(rb[q], gb)
        if pa is None or pb is None:
            continue
        if pa["chrom"] == pb["chrom"] and abs(pa["rs"] - pb["rs"]) < 200:
            continue
        de_tie = abs(pa["de"] - pb["de"]) <= DELTA and max(pa["de"], pb["de"]) <= DE_MAX
        as_tie = min(pa["as_"], pb["as_"]) >= AS_TIE * max(pa["as_"], pb["as_"]) and max(pa["as_"], pb["as_"]) > 0
        n_de += de_tie
        n_as += as_tie
        if de_tie and pa["mapq"] == 0 and pb["mapq"] == 0:
            both_mapq0 += 1
    return dict(shared=len(shared), n_de=n_de, n_as=n_as, both_mapq0=both_mapq0)


def components(loci, edge_fn):
    """connected components over `loci` under boolean edge_fn(i,j). Union-find."""
    parent = {x: x for x in loci}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(len(loci)):
        for j in range(i + 1, len(loci)):
            if edge_fn(loci[i], loci[j]):
                parent[find(loci[i])] = find(loci[j])
    comp = collections.defaultdict(list)
    for x in loci:
        comp[find(x)].append(x)
    return [sorted(v) for v in comp.values()]


def main():
    bam = pysam.AlignmentFile(BAM, "rb")
    # all panel loci as graph vertices; compute pairwise counts once.
    all_loci = sorted({m for _, _, members, _ in PANEL for m in members})
    pc = {}
    for i in range(len(all_loci)):
        for j in range(i + 1, len(all_loci)):
            pc[(all_loci[i], all_loci[j])] = pair_counts(bam, all_loci[i], all_loci[j])

    def key(a, b):
        return (a, b) if (a, b) in pc else (b, a)

    de_edge = lambda a, b: pc[key(a, b)]["n_de"] >= MIN_READS
    as_edge = lambda a, b: pc[key(a, b)]["n_as"] >= MIN_READS

    # OUR DEFINITION: connected components of the de-tie conflict graph.
    de_components = [c for c in components(all_loci, de_edge) if len(c) >= 2]
    as_components = [c for c in components(all_loci, as_edge) if len(c) >= 2]

    # per-candidate ledger: did OUR definition recover it as one component (and not merge it with anything else)?
    print(f"{'family':14} {'regime':34} {'truth':6} {'ours':6} verdict   evidence")
    print("-" * 110)
    rows = []
    tp = tn = fp = fn = 0
    for name, regime, members, truth in PANEL:
        # our definition links members iff they all fall in ONE de-component together (and that component
        # contains exactly these members among the panel — no spurious merge).
        comp_of = {}
        for ci, c in enumerate(de_components):
            for m in c:
                comp_of[m] = ci
        member_comps = {comp_of.get(m) for m in members}
        linked = (len(member_comps) == 1 and None not in member_comps)  # all members in one shared component
        # internal evidence: max de-conflict reads among member pairs
        internal = max((pc[key(a, b)]["n_de"] for a in members for b in members if a < b), default=0)
        verdict = ("TP" if (truth and linked) else "FN" if (truth and not linked)
                   else "FP" if (not truth and linked) else "TN")
        tp += verdict == "TP"; tn += verdict == "TN"; fp += verdict == "FP"; fn += verdict == "FN"
        as_internal = max((pc[key(a, b)]["n_as"] for a in members for b in members if a < b), default=0)
        rows.append(dict(name=name, regime=regime, truth=truth, linked=linked, verdict=verdict,
                         de_reads=internal, as_reads=as_internal))
        print(f"{name:14} {regime:34} {str(truth):6} {str(linked):6} {verdict:9} de={internal} (AS_would={as_internal})")

    print("-" * 110)
    print(f"OUR DEFINITION (de-tie): TP={tp} TN={tn} FP={fp} FN={fn}  "
          f"precision={tp/(tp+fp) if tp+fp else 1.0:.3f} recall={tp/(tp+fn) if tp+fn else 1.0:.3f}")
    print(f"\nFamilies found (de-tie components, >=2 loci): {len(de_components)}")
    for c in de_components:
        print(f"   {{{', '.join(c)}}}")
    print(f"\nFor contrast — AS-tie components: {len(as_components)}")
    for c in as_components:
        spurious = [m for m in c if any(m in mem and not t for _, _, mem, t in PANEL)]
        print(f"   {{{', '.join(c)}}}{'   <- includes AS false-positive loci' if spurious else ''}")

    # explicit FP/FN listing
    print("\n=== FALSE POSITIVES (our definition wrongly links): ===")
    fps = [r for r in rows if r["verdict"] == "FP"]
    print("   none" if not fps else "\n".join(f"   {r['name']} ({r['regime']})" for r in fps))
    print("=== FALSE NEGATIVES (our definition wrongly misses): ===")
    fns = [r for r in rows if r["verdict"] == "FN"]
    print("   none" if not fns else "\n".join(f"   {r['name']} ({r['regime']})" for r in fns))
    print("=== AS-tie false positives our definition AVOIDS: ===")
    avoided = [r for r in rows if not r["truth"] and r["as_reads"] >= MIN_READS and not r["linked"]]
    print("   none" if not avoided else "\n".join(
        f"   {r['name']} ({r['regime']}): AS would link with {r['as_reads']} reads, de correctly rejects ({r['de_reads']})"
        for r in avoided))

    with open(OUT_TSV, "w") as fh:
        fh.write("family\tregime\ttruth_is_family\tour_verdict\tde_conflict_reads\tas_conflict_reads\n")
        for r in rows:
            fh.write(f"{r['name']}\t{r['regime']}\t{int(r['truth'])}\t{r['verdict']}\t{r['de_reads']}\t{r['as_reads']}\n")
    with open(OUT_JSON, "w") as fh:
        json.dump(dict(operating_point=dict(delta=DELTA, de_max=DE_MAX, min_reads=MIN_READS),
                       confusion=dict(tp=tp, tn=tn, fp=fp, fn=fn),
                       de_components=de_components, as_components=as_components, rows=rows,
                       pair_counts={f"{a}~{b}": v for (a, b), v in pc.items() if v["shared"] > 0}), fh, indent=1)
    print(f"\nwrote {OUT_TSV} and {OUT_JSON}")


if __name__ == "__main__":
    main()
