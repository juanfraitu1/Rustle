#!/usr/bin/env python
"""IGV tracks for the copy-assignment result — load these next to the BAM in IGV.

Produces, from a `copy_assign` run + the BAM it swept:
  <out>.tagged.bam (+ .bai) — the reads in the swept regions, each tagged with its assigned copy:
       cp:Z:<family>_c<idx>   (the copy id; IGV: right-click -> "Group by tag" / "Color by tag" -> cp)
       YC:Z:R,G,B             (per-read colour; IGV colours each read by its copy AUTOMATICALLY)
       HP:i:<idx+1>           (so IGV's haplotype colouring also splits the copies)
     A read's SECONDARY placements get the SAME colour, so you SEE the multimapping reads fan out across
     the copies. Tied / ambiguous reads get their own grey, so the K-frontier is visible, not hidden.
  <out>.copies.bed — one line per copy: its genomic extent (from its assigned reads), name = family|copy,
     itemRgb matching the read colour. A loci track to sit above the alignments.

Run:
  python bench/igv_tracks.py --assignments OUT.assignments.tsv --bam reads.bam --regions regions.txt --out OUT
Then in IGV load OUT.tagged.bam and OUT.copies.bed (+ the genome FASTA).
"""
import argparse
import subprocess
from collections import defaultdict

import pysam

# distinct, IGV-friendly colours per copy index; tied/ambiguous get greys (the unresolved classes).
PALETTE = [
    (228, 26, 28), (55, 126, 184), (77, 175, 74), (255, 127, 0), (152, 78, 163),
    (166, 86, 40), (247, 129, 191), (0, 158, 115), (153, 153, 0), (86, 180, 233),
]
TIED_RGB = (150, 150, 150)
AMB_RGB = (200, 200, 200)


def copy_color(idx):
    return PALETTE[idx % len(PALETTE)]


def load_assignments(path):
    """read_name -> (family, copy_idx, status). One row per placement; the call is identical, so first wins."""
    amap = {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[0] in amap:
                continue
            amap[f[0]] = (f[1], int(f[2]), f[3])
    return amap


def tag_for(family, idx, status):
    if status == "assigned":
        return f"{family}_c{idx}", copy_color(idx), idx + 1
    if status == "tied":
        return f"{family}_tied", TIED_RGB, None
    return f"{family}_amb", AMB_RGB, None


def parse_regions(path):
    regs = []
    with open(path) as fh:
        for ln in fh:
            tok = ln.split()
            if not tok:
                continue
            c, rng = tok[0].split(":")
            s, e = rng.split("-")
            regs.append((c, int(s), int(e)))
    return regs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--samtools", default="samtools")
    a = ap.parse_args()

    amap = load_assignments(a.assignments)
    regions = parse_regions(a.regions)
    inbam = pysam.AlignmentFile(a.bam, "rb")
    tmp = f"{a.out}.tagged.unsorted.bam"
    out = pysam.AlignmentFile(tmp, "wb", template=inbam)

    extent = defaultdict(lambda: [None, 1 << 62, 0, None, None])  # key (fam,idx) -> [chrom, start, end, status, idx]
    seen = set()
    n_written = n_tagged = 0
    for (chrom, s, e) in regions:
        for rec in inbam.fetch(chrom, max(0, s), e):
            if rec.is_unmapped:
                continue
            key = (rec.query_name, rec.reference_start, rec.flag)
            if key in seen:
                continue
            seen.add(key)
            info = amap.get(rec.query_name)
            if info:
                fam, idx, status = info
                cp, (r, g, b), hp = tag_for(fam, idx, status)
                rec.set_tag("cp", cp, "Z")
                rec.set_tag("YC", f"{r},{g},{b}", "Z")
                if hp is not None:
                    rec.set_tag("HP", hp, "i")
                n_tagged += 1
                if status == "assigned" and not (rec.is_secondary or rec.is_supplementary):
                    k = (fam, idx)
                    ex = extent[k]
                    ex[0] = rec.reference_name
                    ex[1] = min(ex[1], rec.reference_start)
                    ex[2] = max(ex[2], rec.reference_end or rec.reference_start)
                    ex[3], ex[4] = status, idx
            out.write(rec)
            n_written += 1
    out.close()
    inbam.close()

    sorted_bam = f"{a.out}.tagged.bam"
    subprocess.run([a.samtools, "sort", "-o", sorted_bam, tmp], check=True)
    subprocess.run([a.samtools, "index", sorted_bam], check=True)
    subprocess.run(["rm", "-f", tmp])

    # copies BED (genomic extent of each copy's assigned reads), itemRgb matching the read colour
    with open(f"{a.out}.copies.bed", "w") as fh:
        fh.write('track name="copies" itemRgb="On"\n')
        for (fam, idx), (chrom, st, en, status, ci) in sorted(extent.items()):
            if chrom is None:
                continue
            r, g, b = copy_color(idx)
            fh.write(f"{chrom}\t{st}\t{en}\t{fam}_c{idx}\t0\t+\t{st}\t{en}\t{r},{g},{b}\n")

    print(f"wrote {sorted_bam} (+ .bai): {n_written} alignments, {n_tagged} tagged by copy")
    print(f"wrote {a.out}.copies.bed: {sum(1 for v in extent.values() if v[0])} copies")
    print(f"IGV: load {sorted_bam} + {a.out}.copies.bed + the genome FASTA. Reads auto-colour by copy (YC tag);")
    print(f"     right-click the track -> 'Group by tag' -> cp to stack the copies separately.")


if __name__ == "__main__":
    main()
