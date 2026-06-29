#!/usr/bin/env python3
"""Dump the full minimap2 tag vector for every cross-mapping read on the family panel — the raw material for
evaluating which tags discriminate a GENUINE multimapping conflict (read truly could be from either copy)
from an ARTIFACT (read belongs to one copy; the other placement is spurious cross-alignment).

For each panel pair, for each read with a best-overlap placement on BOTH members via DISTINCT records, emit
one row with every candidate-discriminator tag at the parent (A) and copy (B) placement. One BAM pass; the
analysis agents read the TSV, not the BAM.

Candidate discriminators (per-placement unless noted):
  de    gap-compressed divergence            [adopted winner]
  nmr   NM / query_aligned_len               [adopted corroborator]
  ms    max-scoring-segment score
  s1    best chain score                     s2/s1 = placement uniqueness (parent side)
  s2    2nd-best chain score
  cm    # minimizers on chain
  rl    query length in repetitive seeds
  ts    splice strand (+/-/.)
  mapq  mapping quality                       mapq0-at-both = genuine ambiguity
  nintron  # CIGAR N-gaps (splice structure)  parent-spliced vs copy-unspliced = retrocopy tell
  alnlen   query aligned length
"""
import collections
import pysam

BAM = "/home/juanfra/winloci_scratch/GGO.bam"
OUT = "/mnt/c/Users/jfris/Desktop/Rustle/bench/tag_discriminator_dump.tsv"

GENES = {
    "RABL2B": ("NC_086018.1", 48818440, 48832011), "RABL2A": ("NC_073235.2", 15131653, 15147533),
    "APOBEC3D": ("NC_086018.1", 37023493, 37035944), "APOBEC3F": ("NC_086018.1", 37043278, 37058621),
    "AK6": ("NC_073243.2", 40875017, 40892984), "LOC115934278": ("NC_073227.2", 47415194, 47415981),
    "CCDC196": ("NC_073239.2", 83376335, 83388628), "LOC129526440": ("NC_073238.2", 33024134, 33036450),
    "EEF1A1": ("NC_073229.2", 97608632, 97613071), "LOC109023808": ("NC_073243.2", 97380144, 97381766),
    "CNN2": ("NC_073244.2", 6807132, 6819584), "LOC129524764": ("NC_073237.2", 31580552, 31581609),
    "ADH5": ("NC_073227.2", 119368977, 119386893), "LOC101125574": ("NC_073229.2", 89958139, 89959453),
    "ASDURF": ("NC_073235.2", 92106380, 92111318), "ASNSD1": ("NC_073235.2", 92111059, 92115776),
    "CASP8": ("NC_073235.2", 103850223, 103901051), "FLACC1": ("NC_073235.2", 103898904, 103977183),
    "CREB1": ("NC_073235.2", 110133100, 110220981), "METTL21A": ("NC_073235.2", 110198842, 110241136),
    "GCA": ("NC_073235.2", 64785386, 64835227), "KCNH7": ("NC_073235.2", 64832333, 65309180),
    "GPR39": ("NC_073235.2", 34964263, 35195860), "LYPD1": ("NC_073235.2", 35194087, 35224540),
    "CDPF1": ("NC_086018.1", 44250569, 44262082), "PPARA": ("NC_086018.1", 44159462, 44252981),
    "CCDC188": ("NC_086018.1", 17579900, 17589084), "ZDHHC8": ("NC_086018.1", 17564727, 17581526),
}
PANEL = [
    ("genuine", "RABL2A", "RABL2B"), ("genuine", "AK6", "LOC115934278"), ("genuine", "CCDC196", "LOC129526440"),
    ("artifact_retro", "EEF1A1", "LOC109023808"), ("artifact_retro", "CNN2", "LOC129524764"),
    ("artifact_retro", "ADH5", "LOC101125574"), ("artifact_asym", "APOBEC3D", "APOBEC3F"),
    ("domain", "ASDURF", "ASNSD1"), ("domain", "CASP8", "FLACC1"), ("domain", "CCDC188", "ZDHHC8"),
    ("domain", "CDPF1", "PPARA"), ("domain", "CREB1", "METTL21A"), ("domain", "GCA", "KCNH7"),
    ("domain", "GPR39", "LYPD1"),
]
TAGS = ["de", "ms", "s1", "s2", "cm", "rl", "ts"]


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
        d = dict(chrom=c, rs=r.reference_start, re=r.reference_end,
                 ov=overlap(r.reference_start, r.reference_end, s, e), mapq=r.mapping_quality,
                 nm=r.get_tag("NM") if r.has_tag("NM") else 0, alnlen=qlen,
                 nintron=sum(1 for op, _ in r.cigartuples if op == 3))
        for t in TAGS:
            d[t] = r.get_tag(t) if r.has_tag(t) else ""
        out[r.query_name].append(d)
    return out


def best(recs, gene):
    c, s, e = GENES[gene]
    cand = [x for x in recs if x["chrom"] == c and overlap(x["rs"], x["re"], s, e) > 0]
    return max(cand, key=lambda x: x["ov"]) if cand else None


def main():
    bam = pysam.AlignmentFile(BAM, "rb")
    cols = ["pair", "label", "qname"]
    fields = ["de", "nmr", "ms", "s1", "s2", "cm", "rl", "ts", "mapq", "nintron", "alnlen"]
    for side in ("a", "b"):
        cols += [f"{f}_{side}" for f in fields]
    with open(OUT, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for label, ga, gb in PANEL:
            ra, rb = records(bam, ga), records(bam, gb)
            for q in set(ra) & set(rb):
                pa, pb = best(ra[q], ga), best(rb[q], gb)
                if pa is None or pb is None:
                    continue
                if pa["chrom"] == pb["chrom"] and abs(pa["rs"] - pb["rs"]) < 200:
                    continue  # single straddling record, not a distinct placement

                def row(p):
                    return [p["de"], round(p["nm"] / p["alnlen"], 5), p["ms"], p["s1"], p["s2"],
                            p["cm"], p["rl"], p["ts"], p["mapq"], p["nintron"], p["alnlen"]]
                vals = [f"{ga}~{gb}", label, q] + row(pa) + row(pb)
                fh.write("\t".join(str(v) for v in vals) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
